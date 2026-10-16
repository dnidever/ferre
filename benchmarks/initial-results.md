# Initial executable benchmarks

GNU Fortran 13.3.0; default single precision; included 120-object optical test.
One warmup then three measured runs per configuration. Wall time includes grid
loading, fitting, uncertainties, output and sorting. Same unchanged science code.

| Flags | Threads | Median seconds | Range seconds |
|---|---:|---:|---:|
| -O0 -fopenmp | 1 | 3.423 | 3.324–3.574 |
| -O2 -fopenmp | 1 | 1.870 | 1.869–1.975 |
| -O2 -fopenmp | 2 | 1.268 | 1.167–1.418 |
| -O2 -fopenmp | 4 | 1.218 | 1.218–1.267 |

Optimization gave 1.83x at one thread; four threads with optimization gave
2.81x relative to unoptimized one-thread execution. All runs passed OPFILE
comparisons at atol=rtol=1e-5, including cross-build comparisons. Model spectra
and covariance require additional scientific validation; OPFILE agreement alone
is not a full accuracy assessment.

AMD EPYC 9V74 virtual environment, nine visible CPUs, quota eight CPU equivalents.
Shared-host scheduling and a small optical grid limit generalization to APOGEE.
No GPU was tested. These are CPU feasibility results, not production predictions.
Per-run times, input/executable hashes and column differences are in
initial-results.json. No performance or precision changes were made to Fortran.

Real execution caught two harness/documentation issues: this fork installs
bin/ferre, and normal quit(0) executes STOP 1. Both are handled/documented now.
