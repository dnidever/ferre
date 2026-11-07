Measured November 7, 2026 on an x86_64 Mac, GNU Fortran 11.1.0. Baseline is
FERRE master commit 1120818 (the previous wheel defaults, -O2). Candidate
uses the optimized kernels and -O3 defaults from this branch. Both builds use
OpenMP and -fallow-argument-mismatch, without native CPU flags or fast math.

The commands used `benchmarks/fit_batch.py` with the original regular
APOGEE DR17 M-giant 5x5x5 subset, eight objects, three timed trials, one fitting
thread, and the baseline executable as generator for both builds. The initial
warmup is excluded. Grid and observation hashes match between each pair;
every trial verified recovery of all three fitted labels.

| Fit interpolation | Baseline median | Candidate median | Speedup |
|---|---:|---:|---:|
| Linear (`--inter 1`) | 0.1579 s | 0.1239 s | 1.27x |
| Spline (`--inter 4`) | 15.2590 s | 0.2454 s | 62.18x |

Timing includes startup, loading, fitting, uncertainties, and output. CONT=0,
with no extra resampling or convolution. These spectra are already decompressed,
so this does not measure PCA reconstruction. This small three-label example
is not a benchmark of every HGE workload. Full machine-specific input/output
paths were reduced to executable basenames in these records; binary hashes
remain in the JSON. See `docs/performance.rst` for reproducible commands.
