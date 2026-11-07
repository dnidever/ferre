"""Fit the packaged real APOGEE sparse example, using only the standard library."""
from importlib.resources import files
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    executable = os.environ.get('FERRE_TEST_EXECUTABLE', 'ferre')
    override = os.environ.get('FERRE_SPARSE_GRID_DIR')
    source = Path(override) if override else files('ferre').joinpath('data/apogee-sparse')
    stem = 'apogee-dr17-mgiant-sparse'
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        for suffix in ('.hdr', '.nodes', '.unf', '-weights.json'):
            (work / (stem + suffix)).write_bytes(source.joinpath(stem + suffix).read_bytes())
        nodes = [line for line in (work / (stem + '.nodes')).read_text().splitlines() if not line.startswith('#')]
        assert len(nodes) == 56
        assert (work / (stem + '.unf')).stat().st_size == 56 * 8575 * 4
        control = f"""&LISTA
 NDIM=3
 NOV=0
 NOBJ=1
 SYNTHFILE(1)='{stem}.hdr'
 F_FORMAT=1
 INTER=1
 NTHREADS=1
 CONT=0
 PFILE='start.ipf'
 OFFILE='model.mdl'
/
"""
        def run(text, point):
            (work / 'input.nml').write_text(text)
            (work / 'start.ipf').write_text('star ' + point + '\n')
            result = subprocess.run([executable], cwd=work, capture_output=True, text=True)
            log = result.stdout + result.stderr
            assert result.returncode == 0 or (result.returncode == 1 and 'quit: Run ended' in log), log
            return [float(x) for x in (work / 'model.mdl').read_text().split()]
        truth = [-.62, .8, 3550]
        observation = run(control, ' '.join(map(str, truth)))
        assert len(observation) == 8575 and all(map(math.isfinite, observation))
        assert all(math.isnan(x) for x in run(control, '-.62 .25 3700'))
        weights = json.loads((work / (stem + '-weights.json')).read_text())
        (work / 'weights.dat').write_text(' '.join(map(str, weights)) + '\n')
        (work / 'obs.frd').write_text(' '.join(map(str, observation)) + '\n')
        (work / 'obs.err').write_text(' '.join(['.003'] * 8575) + '\n')
        fit = control.replace('NOV=0', """NOV=3
 INDV=1 2 3
 INIT=0
 ALGOR=1
 STOPCR=1e-9
 SIMP=1e-9
 ERRBAR=1
 FFILE='obs.frd'
 ERFILE='obs.err'
 FILTERFILE='weights.dat'
 OPFILE='fit.opf'""")
        run(fit, '-.75 .7 3520')
        recovered = list(map(float, (work / 'fit.opf').read_text().split()[1:4]))
        assert all(abs(a-b) < tol for a,b,tol in zip(truth,recovered,(.02,.03,10))), recovered
        print('Packaged sparse APOGEE interpolation, exclusion, and fit passed:', recovered)


if __name__ == '__main__':
    main()
