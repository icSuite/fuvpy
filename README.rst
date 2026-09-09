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
``backgroundmodel_BS(..., input_variance="img_variance",
variance_weighted=True)``. This multiplies data-fit rows by
``1/sqrt(img_variance)``. It does not replace the empirical ``dgsigma`` used
for Tukey robust weighting. If ``dampingVal`` is nonzero, its scale must be
reassessed because absolute measurement weights change the balance between
the data and regularization.

``backgroundmodel_BS(..., monotonic=True)`` optionally constrains each
temporal component's spatial B-spline coefficients to be non-decreasing with
the existing ``cos(SZA)/cos(DZA)`` coordinate. The default ``monotonic=False``
keeps the released solver. ``monotonic_solver="normal"`` is the low-cost
operational solver; ``"full"`` retains all weighted image rows as a reference
calculation. A bounded solve is accepted only when SciPy reports
``result.success``; ``max_iter`` remains fixed at 1000. By default a failure
raises. ``monotonic_fallback=True`` instead warns and uses the ordinary solver
for the failed IRLS iteration and all later iterations.

Every BS result includes scalar ``dgsolver_flag``: 0 means the standard solver
was requested, 1 means monotonic solving succeeded throughout, and 2 means a
monotonic failure triggered the explicit standard-solver fallback.

For concise fit monitoring, ``backgroundmodel_BS(..., diagnostics=True)``
adds scalar ``dgfit_iterations`` and ``dgfit_converged`` outputs. They are
absent by default. Complete-orbit WIC comparisons must use an
``img_variance`` saved from ``read_idl(..., measurement_variance=True)``;
variance cannot be reconstructed correctly from an already reflattened WIC
``img`` field. The comparison script verifies the stored ``poisson_offset``
and WIC pre-reflattening provenance attributes before using that variance.


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

    git clone https://github.com/aohma/fuvpy

Use or create a virtual environment with ``conda`` or ``mamba`` where all dependencies are installed.
A working environment is also provided in the repository. Here we create an environment called ``fuvpy-env`` as an example::

    conda env create --name fuvpy-env --file fuvpy/binder/fuvpy-env.yml
    conda activate fuvpy-env

Use pip to install in editable (development) mode::
    
    pip install -e ./fuvpy
