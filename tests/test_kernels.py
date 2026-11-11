"""Compile numerical kernel regressions with the production compiler settings.

Run on a build machine with a Fortran compiler. FERRE's shared global state is
initialized by the test driver, and multiple threads call each optimized kernel
against independent reference results. No benchmark timing threshold is used.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    sources = [root / 'src' / name for name in (
        'share.f90', 'grid_geometry.f90', 'quantized_grid.f90', 'booklib.f90', 'timer.f90',
        'checkstat.f90', 'quit.f90', 'getee.f90', 'spl.f90', 'decompress.f90',
    )]
    sources += [root / 'tests/native/reference_spl.f90', root / 'tests/native/test_kernels.f90']
    compiler = shlex.split(os.environ.get('FC', 'gfortran'))
    flags = shlex.split(os.environ.get('FFLAGS', '-O3 -fopenmp -fallow-argument-mismatch'))
    link_flags = shlex.split(os.environ.get('LDFLAGS', ''))
    with tempfile.TemporaryDirectory(prefix='ferre-kernel-tests-') as temporary:
        executable = Path(temporary) / 'test_kernels'
        subprocess.run(compiler + flags + list(map(str, sources)) + link_flags +
                       ['-o', str(executable)], cwd=temporary, check=True)
        subprocess.run([str(executable)], cwd=temporary, check=True)


if __name__ == '__main__':
    main()
