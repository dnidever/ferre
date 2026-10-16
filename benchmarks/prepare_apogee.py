#!/usr/bin/env python3
"""Download a contiguous 4D DR17 GK PCA slice and prepare deterministic test spectra.
Requires numpy, curl, and a compiled FERRE executable. About 65 MB downloaded.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np

BASE = ('https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/'
        'synspec/marcs/giantisotopes/sgGK_200921nlte_lsfa/')
NAME = 'p_apssgGK_200921nlte_lsfa_012_075'
# Last axis is fastest. Original axes: vmicro, C, N, alpha, M/H, logg, Teff.
SHAPE = (5, 11, 6, 8, 15, 10, 11)
FIXED = (2, 6, 1)  # zero-based indices: log vmicro=.07918, C=0, N=0
START = ((FIXED[0]*11 + FIXED[1])*6 + FIXED[2])*8*15*10*11*900*4
SIZE = 8*15*10*11*900*4
TRUTH = np.array([[.1, -.7, 2.3, 4650], [.25, -1.3, 1.7, 4300],
                  [-.1, .1, 3.1, 5100], [.3, -1.8, 1.1, 4050]])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def slice_header(text):
    """Drop the three fixed leading dimensions in every SYNTH header."""
    blocks = re.findall(r'&SYNTH.*?\n\s*/', text, re.S | re.I)
    if len(blocks) != 4:
        raise ValueError('Unexpected number of SYNTH headers')
    for block in blocks:
        shape = re.search(r'\bN_P\s*=([^\n]+)', block)
        if not shape or tuple(map(int, shape[1].split())) != SHAPE:
            raise ValueError('Unexpected source grid dimensions')
    text = re.sub(r'\bN_OF_DIM\s*=\s*7', 'N_OF_DIM = 4', text)
    for key in ['N_P', 'LLIMITS', 'STEPS']:
        text = re.sub(r'\b' + key + r'\s*=([^\n]+)',
                      lambda m: key + ' = ' + ' '.join(m[1].split()[3:]), text)
    text = re.sub(r'\s*LABEL\((\d+)\)\s*=([^\n]+)',
                  lambda m: '' if int(m[1]) < 4 else
                  '\n LABEL(' + str(int(m[1])-3) + ') =' + m[2], text)
    return text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True, help='New directory outside repository')
    a = p.parse_args()
    exe = a.exe.resolve()
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    original = out / 'source.hdr'
    subprocess.run(['curl', '-fLsS', '--max-time', '180', BASE+NAME+'.hdr',
                    '-o', str(original)], check=True)
    headers = out / 'download-headers.txt'
    binary = out / 'subset.unf'
    subprocess.run(['curl', '-fLsS', '--max-time', '300', '--range',
                    f'{START}-{START+SIZE-1}', BASE+NAME+'.unf',
                    '-D', str(headers), '-o', str(binary)], check=True)
    expected = f'content-range: bytes {START}-{START+SIZE-1}/15681600000'
    if binary.stat().st_size != SIZE or expected not in headers.read_text().lower():
        raise ValueError('Server did not deliver the requested exact binary range')
    (out / 'subset.hdr').write_text(slice_header(original.read_text()))
    (out / 'truth.ipf').write_text('\n'.join('truth%d '%i + ' '.join(map(str,row))
                                           for i,row in enumerate(TRUTH,1))+'\n')
    (out / 'generate.nml').write_text("""&LISTA
 NDIM=4 NOV=0 INTER=3
 SYNTHFILE(1)='subset.hdr' F_FORMAT=1 F_ACCESS=0
 PFILE='truth.ipf' OFFILE='truth.mdl' OPFILE='truth.opf'
 PCACHI=0 PCAPROJECT=0 CONT=0
 /
""")
    with (out / 'generate.log').open('w') as log:
        proc = subprocess.run([str(exe), 'generate.nml'], cwd=out, stdout=log,
                              stderr=subprocess.STDOUT, timeout=300)
    messages = re.findall(r'quit:\s*([^\n]+)', (out/'generate.log').read_text())
    if proc.returncode != 0 and not (proc.returncode == 1 and messages == ['Run ended']):
        raise RuntimeError('FERRE generation failed; inspect generate.log')
    flux = np.loadtxt(out / 'truth.mdl')
    if flux.shape != (4,7514) or not np.isfinite(flux).all() or flux.min() <= 0:
        raise ValueError('Unexpected synthetic spectra')
    case = out / 'case'
    case.mkdir()
    for name in ['subset.hdr','subset.unf']:
        shutil.copy2(out/name,case/name)
    rng = np.random.default_rng(20261016)
    obs, errors, params, truth = [], [], [], []
    for snr in [10,30,100]:
        for i,row in enumerate(flux):
            name = f'snr{snr}_star{i+1}'
            obs.append(row+rng.normal(0,1/snr,row.size))
            errors.append(np.full(row.size,1/snr))
            params.append(name+' 0.0 -0.5 2.5 4500')
            truth.append({'id':name,'snr':snr,'truth':TRUTH[i].tolist()})
    np.savetxt(case/'spectra.frd',obs,fmt='%.8e')
    np.savetxt(case/'spectra.err',errors,fmt='%.8e')
    (case/'spectra.ipf').write_text('\n'.join(params)+'\n')
    (case/'truth.json').write_text(json.dumps(truth,indent=2)+'\n')
    (case/'input.nml').write_text("""&LISTA
 NDIM=4 NOV=4 INDV=1 2 3 4 INTER=3
 SYNTHFILE(1)='subset.hdr' F_FORMAT=1 F_ACCESS=0
 PFILE='spectra.ipf' FFILE='spectra.frd' ERFILE='spectra.err'
 OPFILE='fit.opf' OFFILE='fit.mdl'
 PCACHI=0 PCAPROJECT=0 CONT=0 INIT=0 NRUNS=1 ERRBAR=1 COVPRINT=1
 /
""")
    (case/'provenance.json').write_text(json.dumps({
        'header_url':BASE+NAME+'.hdr','binary_url':BASE+NAME+'.unf',
        'source_header_sha256':sha(original),'slice_sha256':sha(binary),
        'generator_sha256':sha(exe),'byte_start':START,'byte_count':SIZE,
        'original_shape':SHAPE,'fixed_zero_based_indices':FIXED,
        'fixed_physical_values':{'log10_vmicro':.07918,'C':0,'N':0},
        'seed':20261016,'axes':['alpha','M/H','logg','Teff'],
        'limitations':'4D slice, synthetic self-consistency, fixed continuum, not full ASPCAP'
    },indent=2)+'\n')
    print(f'Prepared {case}')


if __name__ == '__main__':
    main()
