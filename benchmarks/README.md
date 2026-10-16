# FERRE benchmarks

This first benchmark measures complete FERRE process wall time and CPU thread
scaling, with numerical regression checks. It does not yet measure isolated
interpolation, model evaluation counts, peak memory, or individual fitting phases.
No scientific Fortran code or production defaults are changed.

## Build and run

Requires Python 3.9+ and a Fortran compiler. Build serially because the existing
makefile does not express all module dependencies. From the repository root:

```bash
make -C src clean
make -C src F90=gfortran OPT='-O3 -march=native -fopenmp'
python3 benchmarks/run.py --exe bin/ferre.x --case test \
  --threads 1 2 4 8 --repeats 3 --out /tmp/ferre-optimized \
  --build-label 'gfortran: record version here; -O3 -march=native -fopenmp'
```

Compile on the target CPU, or replace `-march=native` with a portable target.
Keep the existing precision for initial comparisons. Do not start with fast-math.
To compare an existing production binary against a new build, run it first,
then add `--baseline /path/to/production/results.json` to the candidate command.
Use the same thread list and inputs. Output directories must be new and outside
the input case directory. Large inputs are copied for each run, outside timing;
plan local scratch space accordingly. The original case directory is untouched.

The included test is a small three-dimensional optical-grid smoke/scaling test,
not an APOGEE performance prediction. Its pure-interpolation fixture has no OPFILE
and is intentionally not accepted by this full-fit runner.

## What is recorded

Each thread count has one discarded warmup and three measured runs by default.
JSON includes median/min/max time, objects/hour, executable and input SHA256,
platform, build label, and per-run numerical comparisons. Logs and OPFILE,
OFFILE and SFFILE outputs are retained. The first successful run supplies a
reference unless a baseline is specified. OPFILE rows are aligned by unique
object ID, so output ordering does not affect comparisons. All numerical columns
must be finite and satisfy configurable `--atol`/`--rtol`; failures exit nonzero.
Tolerances apply to the printed columns and are a regression screen, not a
scientific acceptance criterion. Inspect failure/sentinel values, parameter and
abundance shifts, covariance and model spectra before accepting an optimization.
Model spectrum files are saved but not numerically compared automatically.

Whole-process times include startup, grid loading, error estimation, file output,
and sorting. Warmups warm the OS file cache; this is not a cold-disk experiment
or a persistent GPU-grid throughput benchmark. No root cache-flushing is needed.
Threads are requested in both the namelist and environment; check the actual
thread count printed in each log, especially for a binary built without OpenMP.
External BLAS threads are held to one. The runner does not change random seeds:
use deterministic production settings before interpreting differences.

## Production case preparation

Create a dedicated directory containing one `&LISTA` namelist, referenced input
files, and an explicit relative OPFILE (and preferably OFFILE). All output paths
must stay in this directory. Keep large grids/data outside Git. Use an unmodified
production configuration first, including continuum, masks, NRUNS and ERRBAR.
The parser accepts simple quoted paths and single-line comments; use ordinary
relative filenames without exclamation marks. Avoid additional output paths or
scratch settings that point outside the staged directory.

Start with ~20 representative APOGEE spectra: temperature, gravity, metallicity,
S/N and difficult cases. Follow with 100–1000 stars for throughput, on the same
idle compute node. Record CPU model, memory, compiler version, scheduler settings,
filesystem, and build flags alongside results. Repeat on relevant thread counts.

## Next benchmark stages

1. Instrument grid loading, interpolation, PCA reconstruction, continuum,
   optimization and uncertainty estimation separately, plus evaluation counts.
2. Add an interpolation driver with fixed interior/boundary points and both
   nearby trajectories and dispersed points; verify flux agreement with FERRE.
3. Profile a real APOGEE fit before selecting CPU changes or a GPU kernel.
4. Compare GPU startup/transfers and steady-state batches separately. Preserve
   precision and interpolation semantics; measure memory capacity and spectra,
   parameters, uncertainties and failure rates, not just kernel time.

Run harness tests: `python3 -m unittest discover -s benchmarks -p 'test_*.py'`.
