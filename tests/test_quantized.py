"""Native uint16 PCA storage: interpolation, fitting, layout, threads and failures."""

import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import f90nml
from ferre.quantization import quantize_grid, dequantize_grid


class QuantizedGridTests(unittest.TestCase):
    def setUp(self):
        self.executable = os.environ.get("FERRE_TEST_EXECUTABLE") or shutil.which(
            "ferre"
        )
        if not self.executable:
            self.skipTest("Set FERRE_TEST_EXECUTABLE")
        self.executable = str(Path(self.executable).resolve())
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def library(self, sparse=False, transposed=False):
        folder = self.work / "source"
        folder.mkdir()
        axes = (
            ([0, 0.4, 1.3, 2, 4], [0, 0.2, 0.7, 1.2, 2])
            if sparse
            else ([0, 1, 2, 3, 4], [0, 1, 2, 3, 4])
        )
        nodes = list(itertools.product(*axes))
        if sparse:
            nodes = [p for p in nodes if p not in [(0, 2), (4, 0)]]
            nodes = (
                nodes[::2] + nodes[1::2]
            )  # row number must not be inferred from axis products
        self.nodes = nodes
        coefs = np.array(
            [
                [
                    0.1 * x + 0.01 * y,
                    -0.03 * y + 0.01 * x * x,
                    0.125,
                    -0.04 * x + 0.02 * y,
                    0.05 * x * y,
                    0.03 * x * x + 0.02 * y * y,
                ]
                for x, y in nodes
            ],
            dtype="<f4",
        )
        root = dict(
            id="quant-test",
            multi=2,
            n_of_dim=2,
            label=["X", "Y"],
            npix=6,
            npca=[16, 16],
            transposed=int(transposed),
        )
        if sparse:
            root["nodefile"] = "grid.nodes"
            np.savetxt(folder / "grid.nodes", nodes)
        else:
            root.update(n_p=[5, 5], llimits=[0, 0], steps=[1, 1])
        self.basis = np.random.default_rng(17).normal(0, 0.12, (3, 32))
        header = folder / "grid.hdr"
        with header.open("w") as stream:
            stream.write(str(f90nml.Namelist(synth=root)) + "\n")
            for wave in [15000, 16000]:
                stream.write(
                    str(
                        f90nml.Namelist(
                            synth=dict(
                                npix=16,
                                wave=[wave, 1],
                                logw=0,
                                vacuum=1,
                                resolution=22500,
                            )
                        )
                    )
                    + "\n"
                )
            np.savetxt(stream, np.vstack([np.ones(32), np.ones(32), self.basis]))
        (folder / "grid.unf").write_bytes((coefs.T if transposed else coefs).tobytes())
        return header

    def run_grid(
        self, header, access=0, order=1, points=None, threads=1, fit=None, name="model"
    ):
        points = points or [(1.1, 0.8), (1.8, 1.1), (2.2, 0.6)]
        pfile = self.work / (name + ".ipf")
        pfile.write_text(
            "".join(
                f"s{i} " + " ".join(map(str, p)) + "\n" for i, p in enumerate(points)
            )
        )
        config = dict(
            ndim=2,
            nov=0,
            nobj=len(points),
            synthfile=[str(header)],
            f_format=1,
            f_access=access,
            inter=order,
            nthreads=threads,
            cont=0,
            pcachi=0,
            pcaproject=0,
            pfile=str(pfile),
            offile=name + ".mdl",
            opfile=name + ".opf",
        )
        if fit is not None:
            np.savetxt(self.work / "obs.frd", [fit])
            np.savetxt(self.work / "obs.err", [np.full(32, 0.01)])
            config.update(
                nov=2,
                indv=[1, 2],
                init=0,
                algor=1,
                stopcr=1e-9,
                simp=1e-9,
                errbar=0,
                ffile="obs.frd",
                erfile="obs.err",
                opfile=name + ".opf",
            )
        control = self.work / (name + ".nml")
        f90nml.Namelist(lista=config).write(control, force=True)
        result = subprocess.run(
            [self.executable, str(control)],
            cwd=self.work,
            capture_output=True,
            text=True,
        )
        log = result.stdout + result.stderr
        good = result.returncode == 0 or (
            result.returncode == 1 and "quit: Run ended" in log
        )
        if not good:
            raise AssertionError(log)
        model = np.loadtxt(self.work / (name + ".mdl"), ndmin=2)
        return model, log, control

    def test_all_interpolators_access_modes_layouts_and_threads(self):
        for sparse, transposed in itertools.product([False, True], repeat=2):
            with self.subTest(sparse=sparse, transposed=transposed):
                for directory in ["source", "quant", "decoded"]:
                    shutil.rmtree(self.work / directory, ignore_errors=True)
                header = self.library(sparse, transposed)
                quant = quantize_grid(header, self.work / "quant", batch_size=7)
                decoded = dequantize_grid(quant, self.work / "decoded", batch_size=9)
                for order in range(5):
                    ref, _, _ = self.run_grid(decoded, order=order, name="reference")
                    self.assertTrue(np.isfinite(ref[1]).all())
                    for access in (0, 1):
                        actual, log, _ = self.run_grid(
                            quant, access, order, threads=3, name="quant"
                        )
                        np.testing.assert_allclose(
                            actual, ref, atol=1e-6, rtol=0, equal_nan=True
                        )
                        expected = 0 if access else len(self.nodes) * 6 * 2
                        self.assertIn(
                            f"uint16 coefficient array: {expected} RAM bytes", log
                        )
                for directory in ["source", "quant", "decoded"]:
                    shutil.rmtree(self.work / directory)

    def test_fit_and_sparse_exclusion(self):
        header = self.library(sparse=True)
        quant = quantize_grid(header, self.work / "quant")
        decoded = dequantize_grid(quant, self.work / "decoded")
        truth = (1.6, 0.9)
        target, _, _ = self.run_grid(decoded, points=[truth])
        for access in [0, 1]:
            self.run_grid(quant, access, points=[(1.5, 0.8)], fit=target[0], name="fit")
            values = (self.work / "fit.opf").read_text().split()
            np.testing.assert_allclose(
                list(map(float, values[1:3])), truth, atol=0.004, rtol=0
            )
        invalid, _, _ = self.run_grid(quant, points=[(0, 2)], name="excluded")
        self.assertTrue(np.isnan(invalid).all())

    def test_zero_codes_are_valid_pca_mean_spectra(self):
        header = self.library(sparse=True)
        header.with_suffix(".unf").write_bytes(bytes(len(self.nodes) * 6 * 4))
        quant = quantize_grid(header, self.work / "quant")
        for access in (0, 1):
            model, _, _ = self.run_grid(quant, access, order=1)
            np.testing.assert_array_equal(model, np.ones((3, 32)))

    def test_switching_grids_and_access_modes_in_one_process(self):
        header = self.library()
        quant = quantize_grid(header, self.work / "quant")
        decoded = dequantize_grid(quant, self.work / "decoded")
        controls = []
        for i, (path, access) in enumerate(
            [(quant, 0), (quant, 1), (quant, 1), (decoded, 0), (quant, 0)]
        ):
            _, _, control = self.run_grid(path, access, name=f"case{i}")
            controls.append(control)
        listing = self.work / "controls.txt"
        listing.write_text("".join(str(p) + "\n" for p in controls))
        result = subprocess.run(
            [self.executable, "-l", str(listing)],
            cwd=self.work,
            capture_output=True,
            text=True,
        )
        self.assertIn("quit: Run ended", result.stdout + result.stderr)
        self.assertNotIn("already allocated", result.stdout + result.stderr)
        expected = np.loadtxt(self.work / "case0.mdl")
        for i in range(1, 5):
            np.testing.assert_allclose(
                np.loadtxt(self.work / f"case{i}.mdl"), expected, atol=1e-6
            )

    def test_invalid_tables_and_payloads_fail(self):
        header = self.library()
        quant = quantize_grid(header, self.work / "quant")
        table = quant.with_suffix(".quant")
        original = table.read_text()
        for change in ["short", "version", "scale", "missing", "ncoeff"]:
            with self.subTest(change=change):
                table.write_text(original)
                if change == "short":
                    quant.with_suffix(".u16").write_bytes(b"bad")
                elif change == "missing":
                    table.unlink()
                else:
                    n = f90nml.read(table)
                    if change == "version":
                        n["quantization"]["version"] = 99
                    elif change == "ncoeff":
                        n["quantization"]["ncoeff"] = 2
                    else:
                        n["quantization"]["scale"][0] = -1
                    n.write(table, force=True)
                with self.assertRaises(AssertionError) as failure:
                    self.run_grid(quant)
                self.assertIn("FERRE grid error", str(failure.exception))
                if change == "short":
                    shutil.rmtree(self.work / "quant")
                    quant = quantize_grid(header, self.work / "quant")
                    table = quant.with_suffix(".quant")

    def test_converter_does_not_overwrite_and_matches_uint16_roundtrip(self):
        header = self.library()
        quant = quantize_grid(header, self.work / "quant", batch_size=6)
        with self.assertRaises(FileExistsError):
            quantize_grid(header, self.work / "quant")
        decoded = dequantize_grid(quant, self.work / "decoded")
        params = f90nml.read(quant.with_suffix(".quant"))["quantization"]
        q = np.fromfile(quant.with_suffix(".u16"), dtype="<u2").reshape(-1, 6)
        expected = (np.array(params["offset"]) + q * np.array(params["scale"])).astype(
            "<f4"
        )
        np.testing.assert_array_equal(
            np.fromfile(decoded.with_suffix(".unf"), dtype="<f4").reshape(-1, 6),
            expected,
        )
        self.assertEqual(
            quant.with_suffix(".u16").stat().st_size * 2,
            header.with_suffix(".unf").stat().st_size,
        )


if __name__ == "__main__":
    unittest.main()
