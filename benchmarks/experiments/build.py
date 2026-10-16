#!/usr/bin/env python3
"""Build an isolated experimental FERRE executable without editing src/."""
import argparse
from pathlib import Path
import shutil
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--out', type=Path, required=True, help='New build directory outside repository')
p.add_argument('--variant', choices=['original','pca','fused'], default='original')
p.add_argument('--compiler', default='gfortran')
p.add_argument('--flags', default='-O2 -fopenmp')
a = p.parse_args()
root = Path(__file__).resolve().parents[2]
dst = a.out.resolve()
if dst == root or root in dst.parents:
    p.error('Build outside the repository')
dst.mkdir(parents=True, exist_ok=False)
for source in (root/'src').glob('*.f90'):
    shutil.copy2(source,dst/source.name)
shutil.copy2(root/'src/makefile',dst/'makefile')
if a.variant in ['pca','fused']:
    shutil.copy2(Path(__file__).with_name('decompress.f90'),dst/'decompress.f90')
if a.variant == 'fused':
    shutil.copy2(Path(__file__).with_name('cub.f90'),dst/'cub.f90')
subprocess.run(['make','-C',str(dst),'a.out','F90='+a.compiler,'OPT='+a.flags],check=True)
print(dst/'a.out')
