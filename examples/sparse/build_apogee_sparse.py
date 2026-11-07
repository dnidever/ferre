"""Select an uneven, sparse 3D example from the original DR17 5x5x5 subset.

Run download_dr17_mgiant.py first to reconstruct that subset from SDSS, then
pass its output directory and the desired destination to this script.
Only Python's standard library is required for this selection step.
"""
import argparse
from itertools import product
from pathlib import Path
import shutil

STEM = 'apogee-dr17-mgiant-sparse'
HEADER = """&SYNTH
 ID='APOGEE DR17 M giant sparse 3D'
 N_OF_DIM=3
 LABEL(1)='METALS'
 LABEL(2)='LOGG'
 LABEL(3)='TEFF'
 NODEFILE='apogee-dr17-mgiant-sparse.nodes'
 NPIX=8575
 WAVE=4.179000 0.000006
 LOGW=1
 VACUUM=1
 RESOLUTION=22500
/
"""


def build(source, destination):
    """Copy actual source spectra, retaining all 8575 wavelength pixels.

    Select indices 0,1,2,4 on each original five-point axis. Remove the
    cool/high-gravity and hot/low-gravity corner at every metallicity as
    illustrative exclusions, rather than a claim about physical validity.
    Spectrum records and coordinate rows are written in the same order.
    """
    data = (source / 'apogee-dr17-mgiant-full.unf').read_bytes()
    record_bytes = 8575 * 4
    if len(data) != 125 * record_bytes:
        raise ValueError('Expected the original 125-record float32 APOGEE subset')
    destination.mkdir(parents=True, exist_ok=True)
    rows, spectra = [], []
    for m, g, t in product((0, 1, 2, 4), repeat=3):
        if (t, g) in ((0, 4), (4, 0)):
            continue
        record = (m * 5 + g) * 5 + t
        spectra.append(data[record * record_bytes:(record + 1) * record_bytes])
        rows.append(f'{-1 + .25*m:g} {.5*g:g} {3400 + 100*t:g}')
    (destination / f'{STEM}.hdr').write_text(HEADER, encoding='ascii')
    (destination / f'{STEM}.nodes').write_text('# METALS LOGG TEFF\n' + '\n'.join(rows) + '\n', encoding='ascii')
    (destination / f'{STEM}.unf').write_bytes(b''.join(spectra))
    shutil.copyfile(source / 'apogee-dr17-mgiant-full-weights.json', destination / f'{STEM}-weights.json')
    print(f'Wrote {len(rows)} spectra ({len(spectra)*record_bytes:,} bytes)')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    build(args.source, args.destination)
