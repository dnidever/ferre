# Five-dimensional APOGEE nitrogen benchmark

Downloaded **285,120,000 bytes (285.12 MB)** from the public DR17 Synspec GK
LSF-a PCA binary, plus its 14.7 MB header/PCA basis. This is a slice of the
15.7 GB source file, not a full-grid download. Server Content-Range verified:
bytes **7983360000–8268479999 / 15681600000**.

The varying axes are nitrogen (6), alpha abundance (8), metallicity (15),
logg (10), and Teff (11), retaining **79,200 grid spectra**, each represented
by 900 PCA coefficients and reconstructed to 7,514 flux pixels. Carbon is
fixed at zero and log10(vmicro)=0.07918. No Fortran science code was changed.

## Validation

The N=0 plane matches the previous 4D binary slice byte-for-byte. Generating
spectra from the 5D grid at N=0 also reproduces the earlier 4D spectra exactly
at written precision. Six automated harness/header tests pass; the default
4D header generation remains byte-for-byte compatible.

Four off-grid parameter combinations, with N=0.15,0.65,-0.20,0.90, were used
to generate reference spectra. Fixed Gaussian noise at continuum S/N=10,30,100
(seed 20261016) yields 12 synthetic stars. All five parameters are fitted with
cubic Bezier interpolation, ERRBAR=1 and COVPRINT=1. Continuum remains fixed;
INIT=0 uses the same offset guess for every star. One start per fit.

All printed parameters, uncertainties, covariance and other OPFILE columns
matched exactly across builds and thread counts. All fitted model spectra
also matched exactly after alignment by OPFILE object ID. These checks concern
printed output precision and computational reproducibility, not unbiased
recovery or validation of the astrophysical models.

## Timing

GNU Fortran 13.3.0; default single precision. One warmup, three measured runs.
Whole-process wall time includes grid loading, fitting, PCA reconstruction,
uncertainty estimation, output and sorting.

| Compiler flags | Threads | Median seconds | Range seconds |
|---|---:|---:|---:|
| `-O0 -fopenmp` | 1 | 9.893 | 9.542–10.908 |
| `-O2 -fopenmp` | 1 | 2.722 | 2.722–2.973 |
| `-O2 -fopenmp` | 2 | 2.126 | 1.873–2.145 |
| `-O2 -fopenmp` | 4 | 1.470 | 1.419–1.472 |

Optimization alone gave **3.63x** at one thread; optimization with four
threads gave **6.73x** relative to O0/one thread. Runs used the same
shared AMD EPYC 9V74 virtual environment as previous tests. Small batches and
shared-host variability limit the precision of speedup ratios.

## Reproduce

```bash
python3 benchmarks/prepare_apogee.py --dimensions 5 --exe bin/ferre \
  --out /tmp/apogee-5d-fixture
python3 benchmarks/run.py --exe bin/ferre --case /tmp/apogee-5d-fixture/case \
  --nml input.nml --threads 1 2 4 --repeats 3 --out /tmp/apogee-5d-timings
```

The preparation script was run end-to-end against the live archive. Input and
executable hashes, download provenance, and all run times are recorded in
apogee-5d-results.json. Binary assets remain outside Git. Requires numpy and
curl for preparation; the timing runner itself uses the Python standard library.

This adds a fifth interpolation dimension but does not measure full 7D memory
traffic, real-star convergence, production continuum handling, or GPU performance.
