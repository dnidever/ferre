"""Exercise the installed executable, including the APOGEE binary grid."""

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from ferre import get_executable


def run(directory):
    """Run the installed command and check both its status and its output."""
    result = subprocess.run(["ferre"], cwd=directory, capture_output=True, text=True)
    log = result.stdout + result.stderr
    if result.returncode != 0 and not (result.returncode == 1 and "quit: Run ended" in log):
        raise AssertionError(log)
    return log


def main():
    assert shutil.which("ferre")
    assert os.access(get_executable(), os.X_OK)
    assert os.access(get_executable("ascii2bin"), os.X_OK)
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        (work / "grid.dat").write_text("""&SYNTH
 N_OF_DIM=1
 N_P=3
 LABEL(1)='X'
 LLIMITS=0
 STEPS=1
 NPIX=4
 WAVE=15000 1
 LOGW=0
 VACUUM=1
 RESOLUTION=22500
/
1 1 1 1
1 .9 .8 1
1 .8 .6 1
""")
        (work / "start.ipf").write_text("star 0.5\n")
        (work / "input.nml").write_text("""&LISTA
 NDIM=1
 NOV=0
 NOBJ=1
 SYNTHFILE(1)='grid.dat'
 F_FORMAT=0
 INTER=1
 NTHREADS=1
 CONT=0
 PFILE='start.ipf'
 OFFILE='model.mdl'
/
""")
        run(work)
        model = [float(x) for x in (work / "model.mdl").read_text().split()]
        assert len(model) == 4
        assert all(abs(a-b) < 1e-5 for a, b in zip(model, (1, .95, .9, 1))), model
        converted = subprocess.run(
            ["ascii2bin"], cwd=work, input="grid.dat\nunf\n",
            capture_output=True, text=True,
        )
        assert (work / "grid.unf").is_file(), converted.stdout + converted.stderr

        grid_dir = os.environ.get("FERRE_TEST_GRID_DIR")
        if grid_dir:
            source = Path(grid_dir).resolve()
            for suffix in ("hdr", "unf"):
                shutil.copy2(source / f"apogee-dr17-mgiant-full.{suffix}", work / f"apogee.{suffix}")
            weights = json.loads((source / "apogee-dr17-mgiant-full-weights.json").read_text())
            (work / "start.ipf").write_text("star -0.32 1.14 3625\n")
            interpolation = """&LISTA
 NDIM=3
 NOV=0
 NOBJ=1
 SYNTHFILE(1)='apogee.hdr'
 F_FORMAT=1
 INTER=1
 NTHREADS=1
 CONT=0
 PFILE='start.ipf'
 OFFILE='apogee.mdl'
/
"""
            (work / "input.nml").write_text(interpolation)
            run(work)
            observation = [float(x) for x in (work / "apogee.mdl").read_text().split()]
            assert len(observation) == 8575 and all(map(math.isfinite, observation))
            observation[600] += .0001
            observation[4200] -= .0001
            (work / "obs.frd").write_text(" ".join(map(str, observation)) + "\n")
            (work / "obs.err").write_text(" ".join([".003"] * 8575) + "\n")
            (work / "weights.dat").write_text(" ".join(map(str, weights)) + "\n")
            (work / "start.ipf").write_text("star -0.5 1 3600\n")
            fit = interpolation.replace("NOV=0", """NOV=3
 INDV=1 2 3
 INIT=0
 ALGOR=1
 ERRBAR=1
 COVPRINT=1
 FFILE='obs.frd'
 ERFILE='obs.err'
 FILTERFILE='weights.dat'
 OPFILE='fit.opf'""")
            (work / "input.nml").write_text(fit)
            run(work)
            values = (work / "fit.opf").read_text().split()
            metals, logg, teff = map(float, values[1:4])
            assert abs(metals + .32) < .02, values
            assert abs(logg - 1.14) < .03, values
            assert abs(teff - 3625) < 10, values
            print("APOGEE fit", metals, logg, teff)
    print("Installed FERRE interpolation and converter passed")


if __name__ == "__main__":
    main()
