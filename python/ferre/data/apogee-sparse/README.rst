APOGEE DR17 sparse M-giant example
=================================

This is a second example, alongside the original regular 5x5x5 subset in
https://github.com/dnidever/aspcap/tree/master/aspcap/data .
It retains actual spectra from that subset; no new spectra are synthesized.

Labels (coordinate-file column order):

* METALS: -1, -0.75, -0.5, 0 dex
* LOGG: 0, 0.5, 1, 2 dex
* TEFF: 3400, 3500, 3600, 3800 K

At every metallicity, (TEFF, LOGG) = (3400, 2) and (3800, 0) are
omitted. These are illustrative exclusions, not a physical validity mask.
The remaining 56 coordinate rows match the 56 binary spectrum records.
Use INTER=1 for the simplest demonstration near these exclusions.

The .unf file is little-endian float32, 8575 pixels per spectrum, with
spectrum records in the same order as the .nodes file. The wavelength is
10**(4.179 + arange(8575)*0.000006) vacuum Angstroms. All modeled detector
pixels are retained. Unity fills detector gaps and unmodeled edges; use the
accompanying weights (zero there) when fitting.

Source
------

SDSS DR17 SYNspec/MARCS M-giant library:
https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/synspec/marcs/giantisotopes/sgM_200921nlte_lsfa/p_apssgM_200921nlte_lsfa_012_075.hdr
and the matching .unf PCA coefficients. The original subset fixes the first
four source axes at zero-based indices (3, 6, 1, 3): LOG10VDOP, C, N, and
alpha. Its last three axes are METALS, LOGG, and TEFF. PCA spectra were
reconstructed before selecting the records included here.

Reproduce from the repository root (download step requires NumPy)::

    python examples/sparse/download_dr17_mgiant.py /tmp/apogee-full
    python examples/sparse/build_apogee_sparse.py /tmp/apogee-full python/ferre/data/apogee-sparse

Access after pip installation::

    from importlib.resources import files
    grid = files('ferre').joinpath('data/apogee-sparse')
    print(grid.joinpath('apogee-dr17-mgiant-sparse.hdr'))

Copy the header, binary file, and coordinate file together when using them
outside this directory. No download is required after installation.
