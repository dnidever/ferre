"""Measure FERRE kernels with a selected compiler and source tree.

Example: python benchmarks/run.py --flags='-O3 -fopenmp -fallow-argument-mismatch'
Use --source-root to compare another checkout with exactly the same driver.
Timings exclude initialization and file I/O. This is not a whole-fit benchmark.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import shlex
import statistics
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--compiler', default=os.environ.get('FC', 'gfortran'))
    parser.add_argument('--flags', default=os.environ.get('FFLAGS', '-O3 -fopenmp -fallow-argument-mismatch'))
    parser.add_argument('--link-flags', default=os.environ.get('LDFLAGS', ''))
    parser.add_argument('--trials', type=int, default=5)
    parser.add_argument('--scale', type=int, default=1, help='multiply the number of evaluations per trial')
    parser.add_argument('--output', type=Path, default=Path('ferre-benchmark.json'))
    args = parser.parse_args()
    if args.trials < 1 or args.scale < 1:
        parser.error('--trials and --scale must be positive')
    root = args.source_root.resolve()
    sources = [root / 'src' / name for name in (
        'share.f90', 'grid_geometry.f90', 'booklib.f90', 'getee.f90',
        'timer.f90', 'checkstat.f90', 'quit.f90', 'lin.f90', 'qua.f90',
        'cub.f90', 'spl.f90', 'vqbezier.f90', 'vcbezier.f90', 'decompress.f90',
    )]
    sources.append(Path(__file__).resolve().with_name('kernels.f90'))
    compiler = shlex.split(args.compiler)
    results = dict(platform=platform.platform(), machine=platform.machine(),
                   compiler=subprocess.check_output(compiler + ['--version'], text=True).splitlines()[0],
                   flags=args.flags, link_flags=args.link_flags, source_root=str(root),
                   threads=1, trials=args.trials, measurements=[])
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'
    with tempfile.TemporaryDirectory(prefix='ferre-benchmark-') as temporary:
        executable = Path(temporary) / 'benchmark'
        subprocess.run(compiler + shlex.split(args.flags) + list(map(str, sources)) +
                       shlex.split(args.link_flags) + ['-o', str(executable)], cwd=temporary, check=True)
        for method, count, label in (
            (1, 1000, 'linear'), (2, 400, 'quadratic'), (3, 200, 'cubic'),
            (4, 120, 'spline'), (5, 1500, 'pca'),
        ):
            count *= args.scale
            command = [str(executable), str(method), str(count)]
            # Warm the executable and array pages before collecting trials.
            subprocess.check_output(command, cwd=temporary, env=env, text=True)
            timings, checksums = [], []
            for _ in range(args.trials):
                fields = subprocess.check_output(command, cwd=temporary, env=env, text=True).split()
                timings.append(float(fields[2]))
                checksums.append(float(fields[3]))
            median = statistics.median(timings)
            measurement = dict(kernel=label, evaluations=count, seconds=timings,
                               median_seconds=median, microseconds_per_evaluation=median/count*1e6,
                               checksums=checksums)
            results['measurements'].append(measurement)
            print(f'{label:10s} {median/count*1e6:10.1f} us/evaluation')
    args.output.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
