ferre
=====

FERRE is a data analysis code written in FORTRAN90. It matches models to data, taking a
set of observations and identifying the model parameters that best reproduce the data, in a chi-
squared sense. Model predictions are to be given as an array whose values are a function of the
model parameters, i.e. numerically. FERRE holds this array in memory, or in a direct-access
binary file, and interpolates in it. The code returns, in addition to the optimal set of parameters,
their error covariance, and the corresponding model prediction.

The FERRE fortran code was written by Carlos Allende Prieto.  David Nidever repackaged it, added
some Python code to call the Fortran program and added it to Pypi so it can be pip installed.

Installing precompiled FERRE
---------------------------

The Build and test FERRE wheels GitHub Actions workflow produces wheels
for Linux x86-64, macOS Intel, and macOS Apple Silicon. Download and unpack
the artifact for your platform from a successful run, then install its wheel::

    python -m pip install /path/to/ferre-1.0.3-py3-none-PLATFORM.whl

The installed package provides ferre and ascii2bin commands and includes
the required Fortran and OpenMP runtime libraries. Run FERRE in a directory
containing input.nml, or pass a control file as supported by FERRE.

Python callers can locate the bundled executable without importing the
scientific utility module::

    from ferre import get_executable
    print(get_executable())

Building from source
--------------------

Installing from a Git checkout or source archive compiles the executables and
requires GNU Fortran. For a wheel that can be shared with other machines, use
the automated workflow, which bundles runtime libraries and tests the wheel
on a separate runner. Wheels are uploaded as Actions artifacts.

The build accepts FC, FFLAGS, and LDFLAGS. For example::

    python -m pip install build
    python -m build

The existing ferre.ferre scientific module retains its separate scientific
dependencies. The base installation provides the command-line tools and
get_executable using only the Python standard library. Wheels are specific
to their operating system and architecture and support Python 3.9 or later.


Sparse grids with uneven spacing
--------------------------------

FERRE can read an explicit coordinate row for each stored spectrum. It builds
sorted axes and a sparse hash lookup once, then locates interpolation neighbors
using the actual coordinates. Missing combinations and zero spectra exclude
the corresponding interpolation regions.

See docs/sparse-grids.rst for the file format and worked example. From a checkout,
run the example after installing the wheel::

    cd examples/sparse
    ferre

The first row in models.txt is a valid interpolated spectrum. The second row
contains NaN values because a required corner is unavailable.


Performance
-----------

Wheel and source builds default to portable -O3 optimization with OpenMP.
PCA reconstruction traverses contiguous wavelength values, and regular-grid
spline interpolation reuses weights across wavelengths. See
`docs/performance.rst <docs/performance.rst>`_ for build overrides, numerical
checks, and reproducible kernel and complete-fit benchmarks.
