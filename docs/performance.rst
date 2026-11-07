Performance and numerical validation
====================================

Both Python wheel builds and the source makefile use portable optimized GNU
Fortran settings by default::

    -O3 -fopenmp -fallow-argument-mismatch

The wheel build still respects FC, FFLAGS, and LDFLAGS. The makefile respects
command-line F90 and OPT overrides. Build flags do not change the library file
format, optimizer tolerances, interpolation order, or floating-point kind.
The current default real kind is single precision.

To build a debug executable from a source checkout::

    make -C src a.out ascii2bin OPT="-O0 -g -fopenmp -fallow-argument-mismatch -fcheck=all"

To build an optimized wheel::

    python -m pip install build
    python -m build

When changing flags in an existing makefile build, clean its object and module
files first, or use a fresh checkout. Python wheel builds compile in a temporary
directory and do not reuse source-tree objects.

Do not enable -Ofast or -ffast-math for production builds: their assumptions
about finite values conflict with FERRE's NaN and excluded-model checks.
Processor-specific -march=native can be explored for a local build, but is not
used for shared wheels. Those wheels must also run on other processors of the
supported architecture.

Optimized kernels
-----------------

PCA reconstruction traverses each block's wavelength values contiguously inside
the component loop. This matches the layout of the PCA basis and preserves the
component summation order for each pixel. Scaling, the additive constant, and
the mean spectrum are applied as before. This optimization only applies when
FERRE reconstructs PCA spectra; the packaged APOGEE examples are already
reconstructed and therefore do not measure that improvement.

Regular-grid INTER=4 interpolation computes its local spline weights once per
dimension and reuses them at every wavelength. Previously, each wavelength
repeated the same four-point spline calculation. The local abscissae and
boundary shifts are unchanged. Rearranging this linear calculation can change
rounding by a few single-precision ulps. The sparse-grid path already reuses
axis weights and is unchanged by this optimization.

Numerical tests
---------------

On a build machine with GNU Fortran, run::

    python tests/test_kernels.py

The native driver compares optimized PCA reconstruction with a scalar
reference using varying coefficients, component counts, scaling, means, and
unequal block sizes. It compares spline interpolation with the previous
per-pixel implementation for one through five labels, randomized spectra and
coordinates, exact nodes, boundary cells, differing interpolation orders,
and in-memory and direct-access storage. Both kernels are exercised with
three OpenMP threads. Nonfinite outputs and numerical mismatches fail the test;
no timing thresholds are imposed.

The wheel workflow runs these numerical checks on Linux, Intel Mac, and Apple
Silicon, with an additional runtime-checked debug build on Linux. Existing
installation, APOGEE fitting, and sparse-grid integration tests also run against
installed wheels on separate machines.

Reproducible kernel benchmarks
------------------------------

From a checkout with a Fortran compiler, run::

    python benchmarks/run.py --output candidate.json
    python benchmarks/run.py --source-root /path/to/baseline --output baseline.json

The same driver is compiled against each source tree, with the same flags by
default. Override --flags to isolate compiler settings, for example::

    python benchmarks/run.py --flags="-O2 -fopenmp -fallow-argument-mismatch" --output o2.json

Each run records the compiler, platform, settings, individual trial timings,
checksums, and median microseconds per evaluation. Initialization, compilation,
file reads, and output formatting are excluded from the timed kernel sections.
The cases use three labels and 8575 pixels for interpolation, and 12 blocks of
75 coefficients reconstructing 7514 pixels for PCA. Use --scale to increase
trial duration. These are single-thread measurements of individual routines,
not claims about the time for a complete fit.

Complete APOGEE fitting benchmark
---------------------------------

After installing a wheel, benchmark the packaged sparse APOGEE example::

    python benchmarks/fit_batch.py --objects 32 --trials 5 --output sparse-fit.json

The driver generates known spectra, starts fits from offset parameters, and
verifies recovery of all three labels in every trial. Wall-clock timing includes
startup, library loading, model evaluation, uncertainties, and output. An
initial warmup fit is excluded. The defaults use CONT=0 and no extra resampling
or convolution. They represent this small example, not every HGE workload.

For an existing regular APOGEE example and spline interpolation::

    python benchmarks/fit_batch.py --grid-dir /path/to/apogee-grid --stem apogee-dr17-mgiant-full --inter 4 --output regular-spline-fit.json

To compare executables, use a fixed --generator for both runs so the observation
files are identical. The JSON includes hashes of the executable, grid, and
observations as well as all timings. Use --threads for parallel-fit measurements.
The regular example can be regenerated with
examples/sparse/download_dr17_mgiant.py. The packaged sparse example has holes
in its four-coordinate axes and should use INTER=1 for these benchmark points.

Interpretation
--------------

Unoptimized -O0 source builds can be several times slower than optimized builds.
The previous wheels already used -O2, so their incremental gain is smaller and
depends on the workload. PCA and INTER=4 improvements only matter when those
routines contribute to the fit. A small grid's startup and output costs can
obscure interpolation gains. Compare matching inputs and flags on the same
machine, and profile representative scientific fits before generalizing.
