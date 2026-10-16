# FERRE CPU and GPU experiments

These optional variants leave the production `src/` files unchanged. Timings
use the 5D public APOGEE nitrogen fixture: 12 synthetic stars, cubic Bezier,
PCA reconstruction to 7,514 pixels, errors and covariance enabled, fixed
continuum. See ../apogee-5d-results.md for grid provenance and preparation.

## CPU measurements

Fresh O2 reference, GNU Fortran 13.3.0, default single precision; one warmup,
three measured whole-process runs per configuration. Same shared AMD EPYC
9V74 virtual environment, nine visible CPUs / quota eight CPU equivalents.
Grid loading, fit, errors, output and sorting are included.

| Variant | Threads | Median seconds | Strict fit comparison |
|---|---:|---:|---|
| Original O2 | 1 | 3.024 | Reference |
| Original O2 | 4 | 1.371 | Exact printed agreement |
| Cache-friendly PCA, O2 | 1 | 2.423 | Exact printed agreement |
| PCA + gather fusion, O2 | 1 | 2.421 | Exact printed agreement |
| PCA + gather fusion, O2 | 4 | 1.067 | Exact printed agreement |
| Original O3/native | 1 | 1.872 | Fails strict tolerance |
| Original O3/native, contraction off | 1 | 1.871 | Exact printed agreement |
| Original O3/native, contraction off | 4 | 1.220 | Exact printed agreement |
| PCA + gather fusion, O3/native, contraction off | 1 | 0.916 | Exact printed agreement |
| PCA + gather fusion, O3/native, contraction off | 4 | 0.616 | Exact printed agreement |

The combined variant gives **3.30x** relative to the fresh O2 one-thread
reference, or **2.23x** comparing four threads with four threads. Small batch
size, startup costs and shared-host scheduling limit precision/generalization.
The fresh reference differs slightly from earlier runs due to timing variation;
comparisons here use this run's reference rather than the earlier O0 timings.

PCA loop reordering changes traversal order across independent pixels while
preserving accumulation order within each pixel. Fusion reads the first four
neighbor spectra directly into the first interpolation stage, avoiding the
initial full working-array copy for in-memory grids. Its standalone incremental
benefit at O2 is not resolved by these timings. The full workspace declaration
is still present; this experiment reduces copying/touched pages, not the
formal allocation size. Disk-backed interpolation follows the original path.

Plain `-O3 -march=native` produced tiny floating-point changes: maximum fitted
metallicity difference 1e-6 dex and maximum all-column difference 1e-4 (in
covariance), fitted flux difference 1e-5. Disabling contraction restored exact
printed agreement on this fixture. This does not guarantee bitwise equality
on other data, compilers or CPU architectures. No fast-math option was used.

### Profile and validation

An O2 single-thread gprof run attributes about 50% of sampled self time to
vcbezier and 35.6% to PCA decompression; grid reading was 4.8%. Sampling is
approximate and profiling overhead differs from the timed production builds.
The profile records 2,055 objective calls and 2,194 PCA reconstructions.
See profile.txt for full call graph and results.json for hashes, all run times,
per-column differences and validation outcomes.

All PCA, fused, contraction-off and combined fits match the reference's printed
parameters, uncertainties, covariance and fitted model spectra exactly. Separate
74-point checks include random interior points and both near-boundary regions
of every dimension. These variants also match the reference flux output exactly
there. Six automated harness/header tests pass. The 7D grid, real observations,
production continuum/LSF configurations and other interpolation orders are not
validated by this experiment; production defaults remain untouched.

## Build an isolated CPU candidate

From the repository root, after preparing the 5D fixture:

```bash
python3 benchmarks/experiments/build.py --out /tmp/ferre-fast \
  --variant fused --flags='-O3 -march=native -ffp-contract=off -fopenmp'
python3 benchmarks/run.py --exe /tmp/ferre-fast/a.out \
  --case /tmp/apogee-5d-fixture/case --nml input.nml \
  --threads 1 4 --repeats 3 --out /tmp/ferre-fast-results
```

Compile on the intended CPU or use a suitable portable architecture target.
Variants: original, pca, fused (PCA plus interpolation gathering fusion).
Use flags `-O2 -fopenmp -pg` on original to reproduce the profile; run at one
thread and inspect gmon.out with gprof. Timings above used equivalent isolated
builds with the same source variants and compiler flags.

## GPU prototype: prepared, no GPU execution yet

No GPU device nodes, nvidia-smi, nvfortran, or installed CuPy are available in
this environment. **There is no measured GPU speedup.** Installing software
cannot provide missing hardware.

batched_model.py implements cubic interpolation and PCA reconstruction with
NumPy or CuPy. It keeps the grid/basis resident, processes parameter vectors in
batches and uses batched matrix operations for reconstruction. It is a forward
model prototype, not a drop-in optimizer, full fitting pipeline or uncertainty
calculator. Only the prepared 4D/5D float32 PCA fixtures and default reverse-axis
interpolation are supported. No continuum, resampling or additional LSF stage.

The NumPy path was exercised against FERRE for 74 interior/near-boundary points.
Maximum absolute flux difference was **5.12e-6**, passing atol=rtol=1e-5.
Matrix multiplication can change summation order relative to Fortran; exact
bitwise agreement is not expected. **The CuPy path itself remains untested.**

On a CUDA GPU host, install an appropriate CuPy package for its CUDA setup:
https://docs.cupy.dev/en/stable/install.html

```bash
OPENBLAS_NUM_THREADS=1 python3 benchmarks/experiments/validate_model.py \
  --header /tmp/apogee-5d-fixture/case/subset.hdr \
  --exe /tmp/ferre-fast/a.out --backend cupy --batch-size 8 \
  --out /tmp/ferre-gpu-validation
```

Use backend numpy for the CPU validation run. validate_model.py creates fixed
interior/boundary points and a FERRE flux reference automatically. For timings
at different batch sizes, reuse its points/reference with batched_model.py and
--batch-size 1, 8, 32 etc. Use a distinct --out JSON path each time.

Timing reports distinguish grid/context setup, first evaluation (including
kernel compilation), warmed evaluation with resident inputs, and a run including
per-batch input/output transfers. Explicit synchronization brackets GPU timing;
this measures completion, not just launch latency. See CuPy timing guidance:
https://docs.cupy.dev/en/stable/user_guide/performance.html

Next GPU work requires executing and validating on real hardware, measuring
memory and batches, then connecting batched evaluations to FERRE optimizers.
Kernel speed alone cannot establish end-to-end fitting speedup.
