#!/usr/bin/env python3
"""Make reproducible interior/boundary model vectors and compare backend to FERRE."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--header',type=Path,required=True)
p.add_argument('--exe',type=Path,required=True)
p.add_argument('--backend',choices=['numpy','cupy'],default='numpy')
p.add_argument('--batch-size',type=int,default=8)
p.add_argument('--out',type=Path,required=True,help='New directory')
a=p.parse_args()
header=a.header.resolve();exe=a.exe.resolve();out=a.out.resolve()
text=header.read_text();first=re.search(r'&SYNTH.*?\n\s*/',text,re.S|re.I)[0]
def values(key,dtype=float):
    return np.array(re.search(r'\b'+key+r'\s*=([^\n]+)',first)[1].split(),dtype=dtype)
shape=values('N_P',int);low=values('LLIMITS');step=values('STEPS');ndim=len(shape)
if ndim not in (4,5):p.error('Use a prepared 4D/5D fixture')
out.mkdir(parents=True,exist_ok=False)
shutil.copy2(header,out/'subset.hdr')
(out/'subset.unf').symlink_to(header.with_suffix('.unf'))
rng=np.random.default_rng(20261016)
normal=list(rng.uniform(.02,.98,(64,ndim)))
for axis in range(ndim):
    for edge in [1e-5,.99999]:
        vector=np.full(ndim,.4);vector[axis]=edge;normal.append(vector)
points=low+np.array(normal)*step*(shape-1)
np.savetxt(out/'points.txt',points,fmt='%.9g')
(out/'points.ipf').write_text('\n'.join(f'point{i:03d} '+' '.join(f'{x:.9g}' for x in row)
                                       for i,row in enumerate(points))+'\n')
(out/'generate.nml').write_text(f"""&LISTA
 NDIM={ndim} NOV=0 INTER=3
 SYNTHFILE(1)='subset.hdr' F_FORMAT=1 F_ACCESS=0
 PFILE='points.ipf' OFFILE='points.mdl' OPFILE='points.opf'
 PCACHI=0 PCAPROJECT=0 CONT=0
 /
""")
with (out/'generate.log').open('w') as log:
    result=subprocess.run([str(exe),'generate.nml'],cwd=out,stdout=log,
                          stderr=subprocess.STDOUT,timeout=300)
messages=re.findall(r'quit:\s*([^\n]+)',(out/'generate.log').read_text())
if result.returncode and not (result.returncode==1 and messages==['Run ended']):
    raise SystemExit('FERRE reference generation failed')
subprocess.run([sys.executable,str(Path(__file__).with_name('batched_model.py')),
                '--header',str(out/'subset.hdr'),'--points',str(out/'points.txt'),
                '--reference',str(out/'points.mdl'),'--backend',a.backend,
                '--batch-size',str(a.batch_size),'--out',str(out/'model-results.json')],check=True)
