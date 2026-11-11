Native 16-bit PCA Coefficients
==============================

This format stores PCA coefficients as unsigned 16-bit integers, calibrated
separately for every stored coefficient (each wavelength block/component).
The PCA basis, mean spectrum, wavelength segments and node coordinates retain
their existing representations. Existing float32 libraries are unchanged.

In-memory access (F_ACCESS=0) retains an int16 grid, with no full float32 copy.
Direct disk access (F_ACCESS=1) retains only offsets/scales and reads the
requested records. Both decode neighboring coefficients into small temporary
floating-point vectors before calling the existing interpolation routines.
Nearest-neighbor, linear, quadratic, cubic Bezier and cubic spline methods
therefore use the same arithmetic as on a decoded float32 reference library.
Regular and sparse NODEFILE libraries and MULTI wavelength segments are supported.

The coefficient array occupies two rather than four bytes per value in RAM
and on disk. Total process memory also includes coordinates, the sparse lookup,
PCA basis and interpolation/fitting workspaces, so it does not halve exactly.

Create a library
----------------

The converter requires NumPy and f90nml, available through the optional extra::

    python -m pip install 'ferre[quantize]'

From Python:

.. code-block:: python

    from ferre import quantize_grid

    header = quantize_grid("float32-grid/grid.hdr", "quantized-grid")

Or use the command line::

    ferre-quantize quantize float32-grid/grid.hdr quantized-grid

Conversion calibrates ranges over every source record, works in batches, refuses
to overwrite a destination, and leaves the source untouched. It accepts
little-endian float32 binary PCA libraries with a regular grid or NODEFILE.
It preserves coordinate order, wavelength sampling and the original PCA basis.
It adds no new PCA truncation. Integer encoding adds a separate quantization
error that must be validated for the intended abundance precision.

To use an existing observation pipeline, supply the new header and select
F_FORMAT=1 with F_ACCESS=0 or 1. Older FERRE versions do not support the new
header keywords and cannot read the compact payload. For an older executable,
expand a compatible float32 copy in batches::

    ferre-quantize decode quantized-grid/grid.hdr decoded-grid

The Python equivalent is ``ferre.dequantize_grid``. Expansion forfeits the disk
and in-memory storage savings when that float32 copy is used.

Format specification, version 1
-------------------------------

The first SYNTH namelist contains::

    COEFF_FORMAT='uint16'
    QUANTFILE='grid.quant'

COEFF_FORMAT defaults to float32 when omitted. QUANTFILE resolves relative to
the header, or may be absolute. Other SYNTH geometry and PCA metadata are
unchanged. The header ``grid.hdr`` has a companion ``grid.u16`` file; this new
extension prevents old software from mistaking the integers for float32 data.

The payload is a raw stream of little-endian uint16 words with no Fortran
record markers. By default, each spectrum has NPIX consecutive coefficients,
in the same block/component order as its original .unf record. With
TRANSPOSED=1, each coefficient has NTOT consecutive spectrum values instead.
The native reader supports both layouts in memory and on disk; spectrum-major
layout is preferable for efficient direct access. A byte-order conversion is
performed when the host is not little endian.

The parameter sidecar is an ASCII namelist::

    &QUANTIZATION
      VERSION=1
      NCOEFF=1260
      OFFSET=...  ! NCOEFF float64 values
      SCALE=...   ! NCOEFF finite positive float64 values
    /

For each coefficient, ``q = round((value-offset)/scale)`` is an unsigned integer
in 0..65535. Calibration uses the source minimum as offset and
``scale=(maximum-minimum)/65535``. Constant coefficients use scale=1 and q=0.
Reconstruction computes ``offset + scale*q`` in float64 and rounds once to
FERRE's float32 working precision. No saturation or clipping is used by the
converter. Tables with missing/nonfinite values, nonpositive scales, mismatched
counts, unsupported versions, or out-of-range endpoints are rejected. Payload
size must be exactly ``2*NPIX*NTOT`` bytes. A checksum manifest is written by
the converter; the native reader validates layout rather than authenticating files.

This initial implementation quantizes PCA coefficients only, not uncompressed
spectra. It requires one SYNTHFILE; MULTI wavelength segments within it are
supported. FIXFILE is unsupported for these PCA grids, as for existing PCA
libraries. All decoded workspaces are bounded by the interpolation stencil;
no whole-grid decoded cache is built. PCA-space projection/chi-squared options
retain their original meaning; quantization does not establish their accuracy.

Validation
----------

``tests/test_quantized.py`` compares native integer reading with an independently
decoded float32 library. It covers every interpolation option, both access
modes, transposed and normal layouts, shuffled sparse rows, multiple threads,
fit recovery, excluded nodes, repeated control files, and malformed payloads
and calibration tables. It also checks the converter round trip and no-overwrite
behavior. Cross-platform installed-wheel CI runs the same tests.
