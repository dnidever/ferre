"""Integration tests for explicit coordinates, sparse records, and exclusions."""

import itertools
import math
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import tempfile
import unittest


class SparseGridTest(unittest.TestCase):
    def setUp(self):
        self.ferre = os.environ.get("FERRE_TEST_EXECUTABLE") or shutil.which("ferre")
        self.converter = os.environ.get("FERRE_TEST_CONVERTER") or shutil.which("ascii2bin")
        if not self.ferre or not self.converter:
            self.skipTest("Install FERRE or set FERRE_TEST_EXECUTABLE and FERRE_TEST_CONVERTER")
        self.ferre, self.converter = map(lambda p: str(Path(p).resolve()), (self.ferre, self.converter))
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.data = self.work / "data"
        self.data.mkdir()
        self.axes = [[3400, 3500, 3600, 3800, 4000],
                     [-.5, 0, .25, .5, 1], [-1, -.5, 0, .5, 1]]
        self.points = list(itertools.product(*self.axes))

    @staticmethod
    def spectrum(point):
        t, g, m = point
        return [1 + .02 * (t - 3600) / 100, .9 + .1*g, .8 + .08*m, 1] * 8

    def library(self, points=None, zeros=(), single_zero=False, binary=False, regular=False):
        points = self.points if points is None else points
        rows = list(points)
        if not regular:
            random.Random(42).shuffle(rows)
        dimensions = len(rows[0])
        names = ["TEFF", "LOGG", "METALS"][:dimensions]
        header = ("&SYNTH\n ID='sparse-test'\n DATE='2026-11-07'\n"
                  f" N_OF_DIM={dimensions}\n NPIX=32\n WAVE=15000 1\n"
                  " LOGW=0\n VACUUM=1\n RESOLUTION=22500\n")
        if regular:
            header += " N_P=5 5 5\n LLIMITS=3400 -.5 -1\n STEPS=100 .375 .5\n"
        else:
            header += " NODEFILE='grid.nodes'\n"
        for i, name in enumerate(names, 1):
            header += f" LABEL({i})='{name}'\n"
        header += "/\n"
        (self.data / "grid.nodes").write_text("# physical coordinates\n\n" + "".join(
            " ".join(map(str, row)) + " # node\n" for row in rows))
        spectra = [([0]*32 if row in zeros else self.spectrum(row)) for row in rows]
        self.rows = rows
        self.spectra = spectra
        if binary:
            (self.data / "grid.hdr").write_text(header)
            flat = list(itertools.chain.from_iterable(spectra))
            (self.data / "grid.unf").write_bytes(struct.pack(f"={len(flat)}f", *flat))
            return self.data / "grid.hdr"
        text = header
        for row, flux in zip(rows, spectra):
            text += ("0" if single_zero and row in zeros else " ".join(map(str, flux))) + "\n"
        (self.data / "grid.dat").write_text(text)
        return self.data / "grid.dat"

    def control(self, library, points, order=1, access=0, cont=0, threads=1, fit=False):
        (self.work / "start.ipf").write_text("".join(
            f"s{i:03d} " + " ".join(map(str, row)) + "\n" for i, row in enumerate(points)))
        extras = ""
        if fit:
            extras = """ INDV=1 2 3
 INIT=0
 ALGOR=1
 STOPCR=1.e-9
 SIMP=1.e-9
 ERRBAR=1
 COVPRINT=1
 FFILE='obs.frd'
 ERFILE='obs.err'
"""
        text = f"""&LISTA
 NDIM=3
 NOV={3 if fit else 0}
 NOBJ={len(points)}
 NTHREADS={threads}
 SYNTHFILE(1)='{library}'
 F_FORMAT={int(library.suffix==".hdr")}
 F_ACCESS={access}
 INTER={order}
 CONT={cont}
 NCONT=0
 PFILE='start.ipf'
 OPFILE='result.opf'
 OFFILE='model.mdl'
{extras}/
"""
        (self.work / "input.nml").write_text(text)

    def execute(self, succeeds=True):
        result = subprocess.run([self.ferre], cwd=self.work, capture_output=True,
                                text=True, errors="replace", timeout=120)
        log = result.stdout + result.stderr
        if succeeds:
            self.assertTrue(result.returncode == 0 or
                            (result.returncode == 1 and "quit: Run ended" in log), log)
        else:
            self.assertNotIn("quit: Run ended", log)
        return log

    def models(self):
        return [[float(x) for x in row.split()] for row in
                (self.work / "model.mdl").read_text().splitlines()]

    def test_all_orders_nonuniform_and_parallel(self):
        path = self.library()
        queries = [(3550, .1, -.2), (3700, .4, .2), (3900, .8, .7)]
        for order in range(1, 5):
            for threads in (1, 3):
                with self.subTest(order=order, threads=threads):
                    self.control(path, queries, order=order, threads=threads)
                    self.execute()
                    for point, model in zip(queries, self.models()):
                        for actual, expected in zip(model, self.spectrum(point)):
                            self.assertAlmostEqual(actual, expected, delta=2e-5)

    def test_nearest_uses_actual_distances(self):
        path = self.library()
        self.control(path, [(3740, .11, -.2)], order=0)
        self.execute()
        for actual, expected in zip(self.models()[0], self.spectrum((3800,0,0))):
            self.assertAlmostEqual(actual,expected,delta=2e-5)

    def test_missing_corner_and_exact_neighbor(self):
        missing = (3800,.5,.5)
        path = self.library([row for row in self.points if row != missing])
        queries = [(3700,.4,.2),(3600,.25,0),(3800,.5,.5)]
        for cont in (0,1,-1):
            self.control(path,queries,cont=cont)
            self.execute()
            models = self.models()
            self.assertTrue(all(math.isnan(x) for x in models[0]))
            self.assertTrue(all(math.isfinite(x) for x in models[1]))
            self.assertTrue(all(math.isnan(x) for x in models[2]))

    def test_zero_marker_ascii_binary_and_conversion(self):
        missing = (3800,.5,.5)
        for single in (False,True):
            path = self.library(zeros=(missing,),single_zero=single)
            self.control(path,[(3700,.4,.2)])
            self.execute()
            self.assertTrue(all(math.isnan(x) for x in self.models()[0]))
            result = subprocess.run([self.converter],input=str(path)+"\n",
                                    cwd=self.work,capture_output=True,text=True,
                                    errors="replace",timeout=60)
            self.assertTrue(path.with_suffix(".unf").exists(),result.stdout+result.stderr)
            for access in (0,1):
                self.control(path.with_suffix(".hdr"),[(3700,.4,.2),(3500,0,0)],access=access)
                self.execute()
                models = self.models()
                self.assertTrue(all(math.isnan(x) for x in models[0]))
                self.assertTrue(all(math.isfinite(x) for x in models[1]))

    def test_transposed_and_formatted_direct_access(self):
        path = self.library()
        query = (3550,.1,-.2)
        # The normal formatted direct-access format has fixed 12-byte fields.
        flat = "".join(f"{value:12.5E}" for row in self.spectra for value in row)
        path.with_suffix(".fmt").write_text(flat)
        self.control(path,[query],access=1)
        self.execute()
        for actual,expected in zip(self.models()[0],self.spectrum(query)):
            self.assertAlmostEqual(actual,expected,delta=2e-5)
        header=path.read_text().split("/\n",1)[0]+" TRANSPOSED=1\n/\n"
        columns=list(zip(*self.spectra))
        path.write_text(header+"".join(" ".join(map(str,column))+"\n" for column in columns))
        self.control(path,[query])
        self.execute()
        for actual,expected in zip(self.models()[0],self.spectrum(query)):
            self.assertAlmostEqual(actual,expected,delta=2e-5)
        path.with_suffix(".hdr").write_text(header)
        flat=list(itertools.chain.from_iterable(columns))
        path.with_suffix(".unf").write_bytes(struct.pack(f"={len(flat)}f",*flat))
        self.control(path.with_suffix(".hdr"),[query],access=1)
        self.execute()
        for actual,expected in zip(self.models()[0],self.spectrum(query)):
            self.assertAlmostEqual(actual,expected,delta=2e-5)

    def test_pca_zero_coefficients_are_valid(self):
        path = self.library()
        header = path.read_text().split("/\n",1)[0].replace(" NPIX=32"," MULTI=1\n NPCA=32\n NPIX=4")
        header += "/\n&SYNTH\n NPIX=32\n WAVE=15000 1\n LOGW=0\n VACUUM=1\n RESOLUTION=22500\n/\n"
        lines = [" ".join(["1"]*32), " ".join(["0"]*32)]
        lines += [" ".join(str(int(i%4==j)) for i in range(32)) for j in range(4)]
        coeff = [[.02*(t-3600)/100, .1*g, .08*m, 0] for t,g,m in self.rows]
        path.write_text(header+"\n".join(lines)+"\n")
        flat = list(itertools.chain.from_iterable(coeff))
        path.with_suffix(".hdr").write_text(path.read_text())
        path.with_suffix(".unf").write_bytes(struct.pack(f"={len(flat)}f",*flat))
        for access in (0,1):
            self.control(path.with_suffix(".hdr"),[(3600,0,0),(3550,.1,-.2)],access=access)
            self.execute()
            models=self.models()
            self.assertEqual(models[0],[1.0]*32)
            expected=[.99,1.01,.984,1]*8
            for actual,want in zip(models[1],expected):
                self.assertAlmostEqual(actual,want,delta=2e-5)

    def test_multiple_detector_segments(self):
        path=self.library()
        initial=path.read_text().split("/\n",1)[0]
        data=path.read_text().split("/\n",1)[1]
        first=initial.replace(" NPIX=32"," MULTI=2\n NPIX=32")
        segment = "&SYNTH\n NPIX=16\n WAVE={wave} 1\n LOGW=0\n VACUUM=1\n RESOLUTION=22500\n/\n"
        path.write_text(first+"/\n"+segment.format(wave=15000)+segment.format(wave=15100)+data)
        self.control(path,[(3550,.1,-.2)])
        self.execute()
        for actual,expected in zip(self.models()[0],self.spectrum((3550,.1,-.2))):
            self.assertAlmostEqual(actual,expected,delta=2e-5)

    def test_uniform_grid_compatibility(self):
        self.axes=[[3400+100*i for i in range(5)],
                   [-.5+.375*i for i in range(5)],[-1+.5*i for i in range(5)]]
        self.points=list(itertools.product(*self.axes))
        query=(3550,.1,-.2)
        for order in range(5):
            path=self.library(regular=True)
            self.control(path,[query],order=order)
            self.execute()
            legacy=self.models()[0]
            path=self.library()
            self.control(path,[query],order=order)
            self.execute()
            for actual,expected in zip(self.models()[0],legacy):
                self.assertAlmostEqual(actual,expected,delta=2e-5)

    def test_higher_order_support_and_domain_edges(self):
        missing=(3400,0,0)
        path=self.library([row for row in self.points if row!=missing])
        queries=[(3550,.1,-.2),(3500,0,0),(4000,1,1),(4100,1,1)]
        for order in range(1,5):
            self.control(path,queries,order=order)
            self.execute()
            models=self.models()
            self.assertEqual(all(math.isfinite(x) for x in models[0]),order==1)
            self.assertTrue(all(math.isfinite(x) for x in models[1]))
            self.assertTrue(all(math.isfinite(x) for x in models[2]))
            self.assertTrue(all(math.isnan(x) for x in models[3]))

    def test_uncertainty_probe_crossing_hole_is_unavailable(self):
        missing=(3800,.25,0)
        truth=(3600,.25,0)
        path=self.library([row for row in self.points if row!=missing])
        flux=self.spectrum(truth)
        flux[7]+=.0001
        (self.work/"obs.frd").write_text(" ".join(map(str,flux))+"\n")
        (self.work/"obs.err").write_text(" ".join([".003"]*32)+"\n")
        self.control(path,[truth],fit=True)
        self.execute()
        result=(self.work/"result.opf").read_text().split()
        self.assertAlmostEqual(float(result[1]),truth[0],delta=.01)
        self.assertTrue(all(float(x)<0 for x in result[4:7]),result)

    def test_fitting_and_covariance(self):
        path = self.library()
        truth = (3635,.33,-.23)
        flux = self.spectrum(truth)
        flux[7] += .0001
        (self.work / "obs.frd").write_text(" ".join(map(str,flux))+"\n")
        (self.work / "obs.err").write_text(" ".join([".003"]*32)+"\n")
        self.control(path,[(3550,.1,-.1)],fit=True)
        self.execute()
        result = (self.work / "result.opf").read_text().split()
        fitted = list(map(float,result[1:4]))
        for actual,expected,tolerance in zip(fitted,truth,(2,.01,.01)):
            self.assertAlmostEqual(actual,expected,delta=tolerance)
        self.assertTrue(all(float(x)>0 for x in result[4:7]),result)

    def test_invalid_initialization_uses_valid_records(self):
        missing = (3800,.5,.5)
        path = self.library([row for row in self.points if row != missing])
        truth = (3450,-.2,-.7)
        (self.work / "obs.frd").write_text(" ".join(map(str,self.spectrum(truth)))+"\n")
        (self.work / "obs.err").write_text(" ".join([".003"]*32)+"\n")
        self.control(path,[(3700,.4,.2)],fit=True)
        self.execute()
        result = list(map(float,(self.work/"result.opf").read_text().split()[1:4]))
        for actual,expected,tolerance in zip(result,truth,(2,.01,.01)):
            self.assertAlmostEqual(actual,expected,delta=tolerance)

    def test_bad_coordinates_and_record_count(self):
        path = self.library(binary=True)
        original = (self.data/"grid.nodes").read_text()
        for bad in (original + "3400 -.5 -1\n", "3400 0\n", "nan 0 0\n", ""):
            (self.data/"grid.nodes").write_text(bad)
            self.control(path,[(3500,0,0)])
            log = self.execute(succeeds=False)
            self.assertIn("FERRE grid error",log)
        (self.data/"grid.nodes").write_text(original)
        with (self.data/"grid.unf").open("ab") as file:
            file.write(b"1234")
        self.control(path,[(3500,0,0)])
        self.assertIn("record count",self.execute(succeeds=False))


if __name__ == "__main__":
    unittest.main()
