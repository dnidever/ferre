Sparse grids and nonuniform label spacing
=========================================

A library can describe each stored spectrum with an explicit coordinate row.
FERRE derives sorted axis values, physical bounds, and record count from those
rows. Coordinate combinations that are absent are unavailable. This supports
rectilinear axes with uneven spacing and excluded corners or internal holes.

File format
-----------

Add NODEFILE to the library's SYNTH namelist::

    &SYNTH
      N_OF_DIM=2
      LABEL(1)='TEFF'
      LABEL(2)='LOGG'
      NODEFILE='giants.nodes'
      NPIX=32
      WAVE=15000 1
      LOGW=0
      VACUUM=1
      RESOLUTION=22500
    /

The coordinate file contains physical values in LABEL order::

    # TEFF LOGG
    3500 0.5
    3600 0.5
    3500 1.0
    3600 1.0
    3800 1.0
    3600 1.5
    3800 1.5
    3800 2.5

Data row i describes spectrum record i. Comments beginning with # or ! and
blank lines are ignored in the coordinate file. Inline comments are supported.
Coordinates must be finite, with exactly N_OF_DIM numeric values per row.
Duplicate tuples are rejected. Rows may be in any order.

No record-number column or node-count field is required. N_P, LLIMITS, STEPS,
and NTOT in an old header are superseded by the coordinate table. Keep NODEFILE
in the main SYNTH header when using a MULTI library. A relative NODEFILE path
is resolved against the library header's directory, including when the library
is referenced by an absolute path from another working directory.

The spectral data retain FERRE's existing layout. In the normal layout, each
record contains the spectrum described by the matching coordinate row.
ASCII and binary files must contain exactly the corresponding number of values.
Transposed files use the same node order across their columns.

Neighbor lookup
---------------

At load time, FERRE extracts the unique coordinate values on each axis and
constructs a hash lookup from integer axis-index tuples to stored records.
The table contains entries for stored nodes only, with a load factor below 0.5.
The geometry is read-only while fitting threads evaluate models.

For linear interpolation, each axis is searched for the surrounding physical
coordinates. The interval fraction is computed using their actual separation.
The Cartesian product of these intervals identifies the required corners.
FERRE looks up those corners and combines their spectra with multilinear
weights. Corners with exactly zero weight are skipped. Coordinates within a
few single-precision floating-point ulps of a node are snapped to that node.

For the example above, a request at 3650 K and log g 1.25 uses records 4, 5, 6,
and 7. Temperature lies one-quarter of the way from 3600 to 3800 K, and gravity
lies halfway from 1.0 to 1.5. The resulting record weights are 0.375, 0.125,
0.375, and 0.125.

A request at 3550 K and log g 1.25 requires the absent combination 3500 K,
log g 1.5 and is rejected.

Exclusions and interpolation modes
----------------------------------

A required absent record invalidates the evaluation. For uncompressed
libraries, an entirely zero spectrum also marks a record unavailable.
Individual zero pixels remain valid. In spectrum-major ASCII files a line
containing one numeric zero can stand for an entire zero spectrum.
ascii2bin expands that marker while preserving the coordinate reference.

For PCA grids the coordinate table defines availability. Zero coefficient
vectors remain valid because they can describe the mean spectrum. Explicit
coordinate rows should be omitted for excluded PCA nodes.

INTER retains its existing meanings:

* 0 selects the nearest coordinate on each axis, using physical distances.
* 1 performs multilinear interpolation.
* 2 performs quadratic Bezier interpolation using three surrounding coordinates.
* 3 performs cubic Bezier interpolation with four coordinates, retaining the
  quadratic boundary treatment.
* 4 performs local natural cubic-spline interpolation with four coordinates.

Bezier derivative estimates use actual coordinate separations. Sparse uniform
grids agree with the corresponding regular-grid interpolation. Each
nonsingleton axis must have enough coordinates for the selected order.
An exact node request uses only that node on the relevant axis. Otherwise,
higher-order interpolation requires its wider stencil; an unavailable required
node rejects the whole evaluation.

This is a grid of shared coordinate axes. If a coordinate occurs anywhere on
an axis, an absent combination inside another interpolation cell is excluded.
FERRE does not search farther away for replacement corners or triangulate
across absent combinations.

Evaluation, fitting, and uncertainties
--------------------------------------

Invalid model evaluations return before normalization, convolution, and
resampling. Fitting uses a fixed finite penalty in the objective function;
the existing objective's factor-of-50 scaling also applies to this penalty.
Standalone interpolation writes NaN flux for invalid requests. If a fit ends
at an invalid request, parameter/error outputs receive FERRE's existing
-999.999 failure sentinel and the model is NaN.

An invalid initial search location triggers a scan of stored-node free
coordinates with the fixed parameters held at their requested values. A scan
cannot guarantee finding a usable cell if the supplied fixed parameters make
every candidate invalid. Such a fit remains invalid and is flagged at output.

The optimizer continues using affine physical-coordinate normalization.
Bounds and reported parameter/error units come from the true axis endpoints.
Grid-search initialization and ALGOR=0/-1 visit stored node coordinates instead
of a dense Cartesian enumeration. With CHIOUT enabled, sparse diagnostic files
contain the physical node coordinates followed by the log weighted residual
sum for each stored node. Singleton axes can be fixed, but cannot be fitted.

Sparse covariance derivatives use small offsets based on the local axis
spacing. If a derivative probe reaches an excluded region or boundary,
uncertainties are reported unavailable. Zero support is never interpreted as
a precise measurement.

Storage and conversion
----------------------

Sparse coordinates work with memory-loaded ASCII and binary libraries,
formatted direct-access files, unformatted direct-access files, transposed
storage, MULTI detector segments in one library, and PCA libraries.
Direct access to transposed storage reads pixel rows to assemble each requested
spectrum and can be slower than spectrum-major storage.

Use one SYNTHFILE for a sparse library. Detector segments can be combined with
MULTI. Separate SYNTHFILE entries currently require the legacy layout.

Run ascii2bin normally to convert a sparse ASCII library. Keep the coordinate
file with the resulting .hdr/.unf pair; its reference is preserved in the
header. Binary data use the floating-point precision and byte order of the
FERRE build, as in existing FERRE libraries.

Libraries without NODEFILE continue through the legacy interpolation path.
The optional scientific Python utility module's legacy grid-description and
four-label interfaces are separate from this native file-format extension.

Worked example and verification
--------------------------------

The examples/sparse directory includes a complete illustrative two-label
library, coordinate file, input parameter file, and control file. From that
directory, run::

    ferre

The first output spectrum begins with 1.03, 0.925, 0.6, 1.0. The second output
spectrum is invalid because one of its cell corners is absent. These fluxes are
constructed examples for testing.

Run the sparse integration suite against an installed wheel::

    python tests/test_sparse.py -v

To test a source build, set FERRE_TEST_EXECUTABLE and FERRE_TEST_CONVERTER to
its executable paths. The wheel workflow runs both the legacy APOGEE
installation test and this sparse suite on fresh Linux and Mac runners.

Packaged APOGEE example
----------------------

A real three-label APOGEE DR17 M-giant example is installed with the package::

    from importlib.resources import files
    directory = files('ferre').joinpath('data/apogee-sparse')
    header = directory.joinpath('apogee-dr17-mgiant-sparse.hdr')
    print(header)

Set SYNTHFILE(1) to this header path, F_FORMAT=1, NDIM=3, and INTER=1.
Input parameter order is METALS, LOGG, TEFF. For example, ``-0.62 0.8 3550``
is supported, while ``-0.62 0.25 3700`` requires an excluded corner and
returns NaN. A fit uses NOV=3 and INDV=1 2 3. Convert the accompanying
``-weights.json`` array to a whitespace-separated FILTERFILE for fitting;
its zero weights mask detector gaps and unmodeled edge pixels.

The grid has 56 spectra and all 8575 apStar wavelength samples. The shared
axes are METALS=[-1,-0.75,-0.5,0], LOGG=[0,0.5,1,2], and
TEFF=[3400,3500,3600,3800]. At each metallicity the (TEFF, LOGG) corners
(3400,2) and (3800,0) are omitted as illustrative exclusions. The package's
example README records provenance, storage details, and reproduction commands.
The original regular subset remains available in the ASPCAP repository.

``tests/test_apogee_sparse.py`` tests interpolation, rejection of an excluded
cell, and recovery of all three parameters from a generated spectrum. The wheel
workflow runs it against the installed data, outside the source checkout.
