Run `python benchmarks/run.py --output candidate.json` from the repository root.
A Fortran compiler is required; `FC`, `FFLAGS`, and `LDFLAGS` are respected.
Compare with the same command plus `--source-root /path/to/baseline` and
`--output baseline.json`. Use identical flags to isolate algorithm changes, or
change `--flags` to compare compiler settings. Compilation takes place in a
temporary directory; no source-tree object files are reused.

The five single-thread kernel cases measure three-dimensional interpolation
(8575 pixels) and PCA reconstruction (12 blocks, 75 components, 7514 pixels).
The JSON records every timing, checksum, compiler, flags, and platform. Timings
exclude initialization and I/O. These are not complete fitting speedups.

For complete APOGEE fit timings use `python benchmarks/fit_batch.py --help`.
See `docs/performance.rst` for interpretation and validation guidance.
