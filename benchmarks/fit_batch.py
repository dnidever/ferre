"""Time complete three-label APOGEE fits, including startup and output.

The packaged sparse example works by default after installation. For the regular
APOGEE subset use --grid-dir PATH --stem apogee-dr17-mgiant-full. INTER=4 is
particularly useful on that regular grid. Both examples are decompressed, so
this benchmark does not measure the PCA reconstruction optimization.
"""
import argparse
import hashlib
from importlib.resources import files
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile
import time


def executable_path(value):
    path = shutil.which(value)
    if path is None:
        raise ValueError(f'Executable not found: {value}')
    return str(Path(path).resolve())


def run(executable, directory):
    start = time.perf_counter()
    result = subprocess.run([executable], cwd=directory, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    log = result.stdout + result.stderr
    if result.returncode != 0 and not (result.returncode == 1 and 'quit: Run ended' in log):
        raise RuntimeError(log)
    return elapsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default='ferre')
    parser.add_argument('--generator', help='fixed executable for generating identical inputs across comparisons')
    parser.add_argument('--grid-dir', type=Path)
    parser.add_argument('--stem', default='apogee-dr17-mgiant-sparse')
    parser.add_argument('--inter', type=int, choices=(1, 2, 3, 4), default=1)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--objects', type=int, default=32)
    parser.add_argument('--trials', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('ferre-fit-benchmark.json'))
    args = parser.parse_args()
    if min(args.threads, args.objects, args.trials) < 1:
        parser.error('--threads, --objects and --trials must be positive')
    executable = executable_path(args.executable)
    generator = executable_path(args.generator) if args.generator else executable
    source = args.grid_dir.resolve() if args.grid_dir else files('ferre').joinpath('data/apogee-sparse')
    with tempfile.TemporaryDirectory(prefix='ferre-fit-benchmark-') as temporary:
        work = Path(temporary)
        for suffix in ('.hdr', '.unf', '.nodes'):
            entry = source.joinpath(args.stem + suffix)
            if entry.is_file():
                (work / (args.stem + suffix)).write_bytes(entry.read_bytes())
        weights = json.loads(source.joinpath(args.stem + '-weights.json').read_text())
        if len(weights) != 8575:
            raise ValueError('This benchmark expects the full 8575-pixel APOGEE example')
        (work / 'weights.dat').write_text(' '.join(map(str, weights)) + '\n')
        truth = [(-.62 + .02*math.sin(i), .8 + .03*math.cos(i), 3550 + 5*math.sin(i*.7))
                 for i in range(args.objects)]
        (work / 'starts.ipf').write_text(''.join(f'star{i:05d} ' + ' '.join(map(str, p)) + '\n'
                                               for i, p in enumerate(truth)))
        control = f"""&LISTA
 NDIM=3
 NOV=0
 SYNTHFILE(1)='{args.stem}.hdr'
 F_FORMAT=1
 INTER={args.inter}
 NTHREADS=1
 CONT=0
 PFILE='starts.ipf'
 OFFILE='models.mdl'
/
"""
        (work / 'input.nml').write_text(control)
        run(generator, work)
        observations = (work / 'models.mdl').read_text()
        values = [float(x) for x in observations.split()]
        if len(values) != args.objects*8575 or not all(map(math.isfinite, values)):
            raise ValueError('Unsupported interpolation requests for this example grid/order')
        input_hash = hashlib.sha256(observations.encode()).hexdigest()
        (work / 'observations.frd').write_text(observations)
        (work / 'errors.err').write_text((' '.join(['.003']*8575) + '\n')*args.objects)
        (work / 'starts.ipf').write_text(''.join(f'star{i:05d} -.75 .7 3520\n' for i in range(args.objects)))
        control = control.replace('NOV=0', """NOV=3
 INDV=1 2 3
 INIT=0
 ALGOR=1
 STOPCR=1e-8
 SIMP=1e-8
 ERRBAR=1
 FFILE='observations.frd'
 ERFILE='errors.err'
 FILTERFILE='weights.dat'
 OPFILE='fits.opf'""")
        # Generate inputs serially so row order and hashes are independent
        # of the requested number of fitting threads.
        control = control.replace('NTHREADS=1', f'NTHREADS={args.threads}')
        (work / 'input.nml').write_text(control)
        def fit():
            elapsed = run(executable, work)
            rows = [line.split() for line in (work / 'fits.opf').read_text().splitlines() if line.strip()]
            recovered = {row[0]: list(map(float, row[1:4])) for row in rows}
            if len(rows) != args.objects or len(recovered) != args.objects:
                raise AssertionError('Fit did not return one result per object')
            for i, expected in enumerate(truth):
                actual = recovered[f'star{i:05d}']
                if not all(math.isfinite(a) and abs(a-b) < tol
                           for a,b,tol in zip(actual,expected,(.02,.03,10))):
                    raise AssertionError((i, expected, actual))
            return elapsed
        fit()  # Warmup, excluded from the recorded trials.
        timings = [fit() for _ in range(args.trials)]
        median = statistics.median(timings)
        result = dict(platform=platform.platform(), machine=platform.machine(), executable=executable,
                      executable_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
                      generator=generator, grid=args.stem, grid_sha256=hashlib.sha256(
                          (work / (args.stem+'.unf')).read_bytes()).hexdigest(),
                      observations_sha256=input_hash, interpolation=args.inter, threads=args.threads,
                      objects=args.objects, seconds=timings, median_seconds=median,
                      fits_per_second=args.objects/median, verified_recovery=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
        print(f'{args.objects} fits: {median:.3f} seconds; {args.objects/median:.1f} fits/second')


if __name__ == '__main__':
    main()
