Overview
========

The fuvpy package contains code to work with global Far-Ultraviolet (FUV) images of aurora.
At present, the package includes functions to remove background emissions and to detect and model auroral boundaries.

The code uses xarray.Dataset objects to store the images.
Functions to read .idl and .sav files produced by fuview3 and xvis (IMAGE and Polar mission software) into this format are provided.
Other formats need to be converted into a xarray.Dataset with coordinates ['date','row','col'] and data_vars = ['img','mlat','mlon','mlt','glat','glon','dza','sza'].
Here, 'date' is the datetime of the images and 'row' and 'col' refer to the row and column of the camera (2D images).
The data_vars are thus 3D arrays:

- 'img' is the measured intensity of each pixel
- 'mlat', 'mlon' and 'mlt' are the magnetic latitude, longitude and local time of each pixel, respectively
- 'glat' and 'glon' are the geographic latitude and longitude of each pixel, respectively
- 'dza' is the viewing angle of each pixel
- 'sza' is the solar zenith angle of each pixel

If needed, the magnetic coordinates can be calculated from the geographic coordinates (or vice versa) using e.g. apexpy (pip install apexpy).
The solar zenith angle can be calculated based on the subsolar point, which is found using fuvpy.subsol()


Optional detector measurement variance
========================================
``read_idl(..., measurement_variance=True)`` adds ``img_variance`` in the
units of ``img`` squared for WIC, SI12, and SI13. SI variance follows the
Poisson counting relation ``a * img``, where ``a`` is the reconstructed
FUVVIEW3 mission-gain times flat-field multiplier. WIC uses the published
reference compound-Poisson approximation ``10.65 * img``. An explicit
``poisson_offset=0.5`` can be used to avoid zero plug-in variance at very low
counts.

This field describes detector statistics conditional on the supplied
calibration. It excludes calibration uncertainty, atmospheric-response
uncertainty, and background-model error. For WIC, it is calculated before the
optional fuvpy reflattening and retained under the assumption that the
estimated additive reflattening background is deterministic.

The variance can be used explicitly in the BS inversion with
``backgroundmodel_BS(..., variance_weighted=True)``. This multiplies data-fit rows by
``1/sqrt(img_variance)``. It does not replace the empirical ``dgsigma`` used
for Tukey robust weighting. ``dampingVal`` is dimensionless: BS scales the
identity regularization by the median supported positive diagonal of
``GTG`` from the initial weighted design matrix and keep that scale fixed
through IRLS. Historical numerical
damping values therefore need to be retuned.

Directional background model
============================
``backgroundmodel_BS(..., spatial_model="legacy_x_azimuth")`` adds a smooth
directional correction to the ordinary background model. Load IDL files with
``viewing_geometry=True`` so the required solar-view factors are available.
``spencer_taper=True`` limits the directional term using the Spencer
spherical-Chapman profile. The returned dataset stores ``dgmodel_base`` and
``dgmodel_sine`` separately while ``dgmodel`` and ``dgimg`` contain the
combined result. The current production candidate is::

    result = backgroundmodel_BS(
        images,
        spatial_model="legacy_x_azimuth",
        directional_x_knots=[-3.5, 0, 1.5, 3.5],
        temporal_knot_spacing=120,
        directional_temporal_knot_spacing=60,
        spencer_taper=True,
        variance_weighted=True,
        monotonic=True,
        monotonic_fallback=True,
        dampingVal=1e-2
    )

``backgroundmodel_BS(..., monotonic=True)`` optionally constrains each
temporal component's spatial B-spline coefficients to be non-decreasing with
the existing ``cos(SZA)/cos(DZA)`` coordinate. The default ``monotonic=False``
uses the ordinary unconstrained solver. A bounded solve is accepted only when
SciPy reports ``result.success``; ``max_iter`` remains fixed at 1000. With
``monotonic_fallback=False``, failure returns a flagged all-NaN product without
terminating a long processing run. ``monotonic_fallback=True`` instead warns
and uses the ordinary solver for the failed IRLS iteration and all later
iterations.

Every BS result includes scalar ``dgsolver_flag``: 0 means the standard solver
was requested, 1 means monotonic solving succeeded throughout, 2 means a
monotonic failure triggered the ordinary fallback, 3 means monotonic fitting
failed without fallback, and 4 means the ordinary fallback also failed.
``dgfit_iterations`` and ``dgfit_converged`` are always included. Complete-
orbit WIC comparisons must use an
``img_variance`` saved from ``read_idl(..., measurement_variance=True)``;
variance cannot be reconstructed correctly from an already reflattened WIC
``img`` field. The comparison script verifies the stored ``poisson_offset``
and WIC pre-reflattening provenance attributes before using that variance.


Separate 2D temporal and independent backgrounds
================================================

``backgroundmodel_BS()`` retains its existing defaults and one-dimensional
model, including its optional tapered sine asymmetry. The separate
``backgroundmodel_BS2D()`` fits a cubic surface in
``x=cos(SZA)/cos(DZA)`` and the exact stored ``solar_view_sine_factor``.
The latter includes illumination and viewing amplitudes; it is not raw
azimuth or a normalized sine. There is no additional asymmetric term.

Load measurement variance and viewing geometry, then call::

    import fuvpy as fuv

    images = fuv.read_idl(
        filenames, viewing_geometry=True,
        measurement_variance=True, poisson_offset=0.5
    )
    result = fuv.backgroundmodel_BS2D(images, dgsigma_freeze_iteration=2)

The xarray accessor ``images.fuv.backgroundmodel_BS2D(...)`` is also available.
Use ``reflat=False`` when loading SI12/SI13. The new method defaults to
measurement-variance weighting, monotonic x coefficients, dense13-control
x knots, four sine controls, a quadratic120-minute temporal target,
``dampingVal=0.01``, ``tukeyVal=5``, and ``stop=0.01``. ``sKnots`` and
``sine_knots`` configure cubic break points. Sine break points span[-1,1].

Sigma is estimated after every temporal fit through the requested positive
integer iteration, then that per-pixel field freezes. At least one following
weighted solve uses the final estimate and the temporal model converges.
Independent frame fits run automatically, sharing only the final orbit
``sigma²=a*background²+b`` relation. They evaluate it at their own current
backgrounds and start with fresh Tukey weights. No temporal shape prior or
inherited pixel weights constrain the frame fits.

Temporal results retain ``dgmodel``, ``dgimg``, ``dgweight`` and ``dgsigma``.
Independent results use ``dgmodel_independent``, ``dgimg_independent``,
``dgweight_independent`` and ``dgsigma_independent``. Both coefficient tensors,
sigma-estimation history, iteration counts and convergence diagnostics are
saved. Temporal fields are never overwritten by the independent stage.
Sigma is empirical residual spread, separate from detector measurement
variance and fitted-background uncertainty. A failed/nonconverged2D solve
raises an error rather than silently returning partial fits.

Shared metadata pointing correction
===================================

Both methods accept a sensor-keyed mapping and optional
``pointing_correction=True``. Load all three sensors with
``pointing_metadata=True``; this retains original spacecraft/camera metadata,
native coordinate grids and millisecond timestamps, and enables viewing
geometry. Existing loader defaults and outputs are unchanged when disabled::

    sensors = {
        'WIC': fuv.read_idl(wic_files, pointing_metadata=True,
                            measurement_variance=True, poisson_offset=0.5),
        'SI12': fuv.read_idl(si12_files, pointing_metadata=True, reflat=False,
                             measurement_variance=True, poisson_offset=0.5),
        'SI13': fuv.read_idl(si13_files, pointing_metadata=True, reflat=False,
                             measurement_variance=True, poisson_offset=0.5)
    }
    results = fuv.backgroundmodel_BS2D(
        sensors, pointing_correction=True,
        pointing_options={'method': 'clock_phase'},
        dgsigma_freeze_iteration=2
    )

``fuv.backgroundmodel_BS(sensors, pointing_correction=True, ...)`` enables the
same correction for the1D method, retaining that method's temporal-only
behavior. A single-sensor correction request is rejected: the shared state
must be repaired consistently in WIC, SI12 and SI13. All corrections are off
by default.

The shared WIC reference repairs isolated spacecraft-position outliers using
a robust local cubic. ``method='clock_phase'`` additionally smooths the
clock-aware phase with an eleven-frame local quadratic; ``method='position'``
selects position-only repair. Camera alignments remain fixed. Corrections
are evaluated at each camera's precise time on its native grid, without
brightness or limb fitting. Geographic coordinates, SZA/DZA, both solar-view
factors and historical magnetic coordinates are recalculated together. The
packaged epoch1997 Apex table preserves FUVVIEW3's coordinate convention.

Counts and detector variance remain unchanged. Corrected science coordinates
retain the intersection of original mapping and corrected visibility/DZA
support. Diagnostics record original spacecraft state, coordinate hashes,
position/phase corrections and unsupported intervals. Cameras associate by
nearest millisecond observation time rather than equal SPIN counters. WIC
counter discontinuities and time gaps separate smoothing segments; short
segments and reference-time boundaries are not extrapolated. This repairs
metadata consistency; it is not an independent attitude determination.

For corpus processing, ``read_idl(..., pointing_metadata=True,
skip_unmapped=True)`` excludes files with neither stored nor reconstructed
valid DZA pixels. Each excluded filename is warned and retained in the
``skipped_unmapped_files`` attribute. Other geometry validation errors still
raise. The default remains strict and does not skip files.

Use complete-orbit metadata before selecting science frames. The public
``fuv.correct_pointing(sensors, ...)`` returns corrected copies separately
when hemisphere selection or other preprocessing must occur before fitting.
After that shared step, call either background method with correction off.
Already corrected inputs reject a second correction request.


Dependencies
============
fuvpy depends on the following:

- matplotlib
- netcdf4
- numpy
- pandas
- scipy
- xarray
- polplot (github/klaundal/polplot)

You should also have git version >= 2.13

Install
=======
Clone the repository::

    git clone https://github.com/icSuite/fuvpy.git

Use an existing scientific Python environment and install in editable mode::
    
    python -m pip install -e ./fuvpy
