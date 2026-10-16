# APOGEE DR17 GK slice benchmark

Actual public data from [sgGK_200921nlte_lsfa](https://data.sdss.org/sas/dr17/apogee/spectro/speclib/synth/synspec/marcs/giantisotopes/sgGK_200921nlte_lsfa/).
The 15,681,600,000-byte FERRE binary is paired with a 14,748,879-byte header.
HTTP range extraction downloaded exactly bytes 8030880000–8078399999 (47,520,000
bytes); the server returned the expected Content-Range and total size.

The source is 7D, with 900 PCA coefficients and 7,514 reconstructed spectral
pixels. We took a contiguous 4D slice: alpha abundance (8), metallicity (15),
logg (10), Teff (11). Microturbulence is fixed at log10(vmicro)=0.07918,
C=0 and N=0. Header axes are reduced consistently in all four SYNTH blocks;
the PCA basis and reconstructed wavelength sampling are preserved.

Four off-grid truth spectra were generated with reference FERRE, each with
Gaussian continuum noise at S/N=10,30,100 (seed 20261016), giving 12 fits. Cubic
Bezier interpolation, fixed normalized continuum, four fitted parameters,
INIT=0 at the same deliberately offset initial guess, one start, ERRBAR=1,
COVPRINT=1. This tests self-consistency and computational regression, not
astrophysical accuracy or the full ASPCAP workflow.

| GNU Fortran 13.3.0 flags | Threads | Median seconds | Range seconds |
|---|---:|---:|---:|
| -O0 -fopenmp | 1 | 5.332 | 4.930–5.389 |
| -O2 -fopenmp | 1 | 1.268 | 1.268–1.319 |
| -O2 -fopenmp | 2 | 0.817 | 0.817–0.967 |
| -O2 -fopenmp | 4 | 0.717 | 0.667–0.718 |

One warmup and three measured runs. Whole-process wall times include grid
loading, PCA reconstruction, fitting, uncertainty calculation, output and
sorting. Default single precision; no science code changes. Same AMD EPYC
9V74 shared virtual environment as initial-results.md. Optimization alone is
4.2x faster; optimization plus four threads is 7.4x relative to O0/one thread.
Shared-host timing and small batch size limit precision of these ratios.

All printed OPFILE numeric columns (including fitted parameters, uncertainties
and covariance) matched exactly across all runs/builds/threads. Fitted model
spectra also matched exactly after alignment by OPFILE object ID. This does not
establish agreement at more precision than FERRE writes to its output files.

## Reproduce

Python requires numpy; curl and the compiled FERRE binary must be available.
From the repository root:

```bash
python3 benchmarks/prepare_apogee.py --exe bin/ferre --out /tmp/apogee-fixture
python3 benchmarks/run.py --exe bin/ferre --case /tmp/apogee-fixture/case \
  --nml input.nml --threads 1 2 4 --repeats 3 --out /tmp/apogee-timings
```

The preparation script was exercised with the verified downloads cached locally
and reproduced the benchmark header, binary, spectra, errors, initial guesses
and namelist byte-for-byte. Hashes and extraction provenance are recorded.
Large assets stay outside Git. Results above used the identical fixture before
adding provenance metadata. The download requests were verified against the
live archive, including the binary range response.

## Limits and next steps

This 4D slice cannot predict the full 7D cubic interpolation workspace cost.
It is suitable for PCA/flux reconstruction, numerical regression and an initial
APOGEE benchmark. Next: the full 7D grid (disk access or a larger-memory host),
real observed spectra, production continuum/masks, larger batches, and timing
instrumentation before a GPU implementation. No GPU benchmark has run.
