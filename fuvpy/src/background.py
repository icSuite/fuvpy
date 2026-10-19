#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:28:25 2022

@author: aohma
"""
from pathlib import Path
from collections.abc import Mapping
import warnings

import numpy as np
import xarray as xr

from scipy.interpolate import BSpline
from scipy.optimize import curve_fit, lsq_linear
from scipy.stats import binned_statistic
from scipy.linalg import lstsq

from fuvpy.utils.sunlight import subsol


def _monotonic_bs_transform(n_scp, n_tcp):
    """Map free baselines and nonnegative increments to BS coefficients.

    The existing coefficient vector is ordered with spatial basis outermost
    and temporal basis innermost. Each temporal component receives one free
    first spatial coefficient followed by cumulative spatial increments.
    """

    n_coefficients = n_scp * n_tcp
    transform = np.zeros((n_coefficients, n_coefficients))

    for temporal in range(n_tcp):
        for spatial in range(n_scp):
            row = spatial*n_tcp + temporal
            for increment in range(spatial + 1):
                column = increment*n_tcp + temporal
                transform[row, column] = 1

    return transform


def _monotonic_base_transform(n_base_scp, n_base_tcp, n_directional_coefficients):
    """Constrain only the leading one-dimensional base coefficients.

    Coefficients after the base belong to the signed directional model and
    remain free.
    """

    n_base_coefficients = n_base_scp*n_base_tcp
    n_coefficients = n_base_coefficients + n_directional_coefficients
    transform = np.eye(n_coefficients)
    lower = np.full(n_coefficients, -np.inf)

    for temporal in range(n_base_tcp):
        baseline = temporal
        for spatial in range(n_base_scp):
            row = spatial*n_base_tcp + temporal
            transform[row, :] = 0
            transform[row, baseline] = 1

            for increment in range(1, spatial + 1):
                column = increment*n_base_tcp + temporal
                transform[row, column] = 1
                lower[column] = 0

    return transform, lower


def _solve_monotonic_bs(weighted_design, weighted_data, damping, n_scp, n_tcp, solver, transform=None, lower=None):
    """Solve a B-spline fit using bounded transformed coefficients.

    BVLS is allowed at most 1000 iterations. SciPy's ``result.success`` is the
    convergence criterion. Without an explicit transform and bounds, every
    spatial coefficient sequence is constrained to be non-decreasing.
    """

    if transform is None:
        transform = _monotonic_bs_transform(n_scp, n_tcp)
    n_coefficients = transform.shape[0]

    if solver == 'full':
        # Reference calculation against every weighted image row.
        design = weighted_design @ transform
        data = weighted_data

        damping_diagonal = np.diag(damping)
        if np.any(damping_diagonal > 0):
            regularization = np.sqrt(damping_diagonal)[:, None] * transform
            design = np.vstack([design, regularization])
            data = np.r_[data, np.zeros(n_coefficients)]
    elif solver in ['normal', 'normal_system']:
        # Operational calculation: reduce all pixels to the small normal
        # system before applying the same coefficient bounds.
        if solver == 'normal':
            normal_matrix = weighted_design.T@weighted_design + damping
            normal_data = weighted_design.T@weighted_data
        else:
            normal_matrix = weighted_design + damping
            normal_data = weighted_data

        normal_matrix = transform.T@normal_matrix@transform
        normal_data = transform.T@normal_data

        try:
            lower_cholesky = np.linalg.cholesky(normal_matrix)
            design = lower_cholesky.T
            data = np.linalg.solve(lower_cholesky, normal_data)
        except np.linalg.LinAlgError:
            symmetric_matrix = (normal_matrix + normal_matrix.T)/2
            eigenvalue, eigenvector = np.linalg.eigh(symmetric_matrix)
            tolerance = (np.max(np.abs(eigenvalue)) * n_coefficients * np.finfo(float).eps)
            if np.min(eigenvalue) < -tolerance:
                raise RuntimeError("Monotonic normal matrix is not positive semidefinite")

            eigenvalue = np.clip(eigenvalue, 0, None)
            design = np.sqrt(eigenvalue)[:, None] * eigenvector.T
            data = np.linalg.lstsq(design.T, normal_data, rcond=None)[0]
    else:
        raise ValueError("solver must be 'full', 'normal', or 'normal_system'")

    if lower is None:
        lower = np.zeros(n_coefficients)
        lower[:n_tcp] = -np.inf

    try:
        result = lsq_linear(design, data, bounds=(lower, np.inf), method='bvls', max_iter=1000)
    except Exception as error:
        raise RuntimeError(f"Monotonic B-spline solver failed: {error}") from error
    if not result.success:
        raise RuntimeError(f"Monotonic B-spline fit failed (optimality={result.optimality:.3g}): {result.message}")

    return transform @ result.x


def _normal_matrix_scale(normal_matrix):
    """Return the median scale of coefficient columns supported by the data."""

    diagonal = np.diag(normal_matrix)
    largest = np.max(diagonal)
    tolerance = largest*diagonal.size*np.finfo(float).eps
    supported = diagonal[diagonal > tolerance]
    if supported.size == 0:
        raise ValueError("The weighted design matrix contains no supported coefficients")

    return np.median(supported)


def _temporal_spline_basis(time, knot_spacing, order):
    """Build a temporal B-spline basis from a target spacing in minutes."""

    if isinstance(knot_spacing, (bool, np.bool_)) or not np.isscalar(knot_spacing) or not np.isfinite(knot_spacing) or knot_spacing <= 0:
        raise ValueError("temporal knot spacing must be a finite positive number of minutes")
    if not isinstance(order, (int, np.integer)) or order < 0:
        raise ValueError("tOrder must be a non-negative integer")

    time = np.asarray(time, dtype=float)
    duration = time[-1] - time[0]
    if duration < 0 or np.any(np.diff(time) < 0):
        raise ValueError("image dates must be ordered from earliest to latest")

    if duration == 0:
        breaks = np.array([time[0], time[0] + 1])
        knots = breaks
        effective_order = 0
        actual_spacing = 0.0
    else:
        n_intervals = max(1, int(np.ceil(duration/knot_spacing)))
        breaks = np.linspace(time[0], time[-1], n_intervals + 1)
        knots = np.r_[np.repeat(breaks[0], order), breaks, np.repeat(breaks[-1], order)]
        effective_order = order
        actual_spacing = duration/n_intervals

    n_control_points = len(knots)-effective_order-1
    design = BSpline(knots, np.eye(n_control_points), effective_order)(time)
    return breaks, design, n_control_points, effective_order, actual_spacing


def _accumulate_normal_system(spatial_designs, frame_support, base_temporal_design, directional_temporal_design, data, robust_weight, spatial_weight, measurement_weight, pixels_per_frame):
    """Accumulate the weighted normal system without a full design matrix.

    Base coefficients come first, followed by directional coefficients. Within
    each block, the spatial basis is outermost and the temporal basis is
    innermost. The full design is assembled one frame at a time, avoiding a
    whole-orbit pixels-by-coefficients matrix.
    """

    n_base_coefficients = spatial_designs[0][0].shape[1]*base_temporal_design.shape[1]
    n_directional_coefficients = 0
    if directional_temporal_design is not None:
        n_directional_coefficients = spatial_designs[0][1].shape[1]*directional_temporal_design.shape[1]
    n_coefficients = n_base_coefficients + n_directional_coefficients
    normal_matrix = np.zeros((n_coefficients, n_coefficients))
    normal_data = np.zeros(n_coefficients)

    for frame, ((base_spatial, directional_spatial), pixel_index) in enumerate(zip(spatial_designs, frame_support)):
        if pixel_index.size == 0:
            continue

        data_index = frame*pixels_per_frame + pixel_index
        fit_weight = np.sqrt(robust_weight[data_index]) * spatial_weight[data_index] * measurement_weight[data_index]
        base_design = (base_spatial[:, :, None]*base_temporal_design[frame][None, None, :]).reshape((pixel_index.size, -1))
        if directional_temporal_design is not None:
            directional_design = (directional_spatial[:, :, None]*directional_temporal_design[frame][None, None, :]).reshape((pixel_index.size, -1))
            frame_design = np.hstack([base_design, directional_design])
        else:
            frame_design = base_design

        weighted_design = frame_design*fit_weight[:, None]
        weighted_data = data[data_index]*fit_weight
        normal_matrix += weighted_design.T@weighted_design
        normal_data += weighted_design.T@weighted_data

    return normal_matrix, normal_data


def _spencer_taper_weights(sza):
    """Return the normalized spherical-Chapman profile used by Lompe.

    The taper is interpolated in SZA and set to zero beyond the table's
    120-degree limit. It limits directional freedom; it is not an IMAGE
    brightness calibration.
    """

    sza = np.asarray(sza, dtype=float)
    table = np.loadtxt(Path(__file__).parents[1]/'data'/'spencer_chapman_taper.txt')
    weight = np.interp(sza, table[:, 0], table[:, 1], left=1, right=0)
    return np.where(np.isfinite(sza), weight, np.nan)


def _nonaxisymmetric_sh_keys(max_degree, max_order):
    """Return the non-axisymmetric terms used by the historical SH model.

    The old residual correction retained terms for which ``n-m`` is even. Its
    cosine expansion also included ``m=0``. Here both cosine and sine parts
    start at ``m=1`` because an axisymmetric SH contribution would compete
    directly with the coestimated one-dimensional base model.
    """

    if not isinstance(max_degree, (int, np.integer)) or max_degree < 1:
        raise ValueError("sh_degree must be a positive integer")
    if not isinstance(max_order, (int, np.integer)) or max_order < 1:
        raise ValueError("sh_order must be a positive integer")

    return [(degree, order) for degree in range(1, max_degree + 1) for order in range(1, min(degree, max_order) + 1) if (degree - order) % 2 == 0]


def _schmidt_legendre(colatitude, max_degree, max_order, keys):
    """Evaluate the historical Schmidt semi-normalized Legendre functions."""

    theta = np.deg2rad(np.asarray(colatitude, dtype=float).reshape(-1, 1))
    sine = np.sin(theta)
    cosine = np.cos(theta)

    functions = {(degree, order): np.zeros_like(theta) for degree in range(max_degree + 1) for order in range(max_degree + 1)}
    functions[(0, 0)] = np.ones_like(theta)
    normalization = {(0, 0): 1.0}

    for degree in range(1, max_degree + 1):
        for order in range(min(degree, max_order) + 1):
            if degree == order:
                value = sine*functions[(degree - 1, order - 1)]
            elif degree == 1:
                value = cosine*functions[(degree - 1, order)]
            else:
                coefficient = ((degree - 1)**2 - order**2)/((2*degree - 1)*(2*degree - 3))
                value = cosine*functions[(degree - 1, order)] - coefficient*functions[(degree - 2, order)]

            if order == 0:
                normalization[(degree, 0)] = normalization[(degree - 1, 0)]*(2*degree - 1)/degree
            else:
                normalization[(degree, order)] = normalization[(degree, order - 1)]*np.sqrt((degree - order + 1)*(2 if order == 1 else 1)/(degree + order))

            functions[(degree, order)] = value

    return np.hstack([functions[key]*normalization[key] for key in keys])


def _real_sh_design(glat, solar_relative_longitude, max_degree, max_order):
    """Build the non-axisymmetric real SH basis used for the joint model."""

    keys = _nonaxisymmetric_sh_keys(max_degree, max_order)
    legendre = _schmidt_legendre(90 - np.asarray(glat), max_degree, max_order, keys)
    longitude = np.deg2rad(np.asarray(solar_relative_longitude)).reshape(-1, 1)
    orders = np.asarray([order for _, order in keys])[None, :]

    cosine = legendre*np.cos(longitude*orders)
    sine = legendre*np.sin(longitude*orders)
    return np.hstack([cosine, sine]), keys


def backgroundmodel_BS(imgs,**kwargs):
    '''
    Function to model the FUV dayglow and subtract it from the input image

    Parameters
    ----------
    imgs : xarray.Dataset or sensor-keyed mapping of Datasets
        FUV images imported by read_idl(). Shared pointing correction requires
        all three sensors in one mapping.
    sOrder : int, optional
        Order of the spatial spline fit. The default is 3.
    dampingVal : float, optional
        Dimensionless damping (zeroth-order Tikhonov regularization). The
        identity matrix is scaled by `dampingVal` times the median supported
        positive diagonal of the initial weighted `GTG`. The default is 0.
    tukeyVal : float, optional
        Determines to what degree outliers are down-weighted.
        Iterative reweights is (1-(residuals/(tukeyVal*rmse))^2)^2
        Larger tukeyVal means less down-weight
        Default is 5
    stop : float, optional
        When to stop the iteration. The default is 0.001.
    sKnots : array like, optional
        Location of the spatial Bspline knots. The default is None (default is used).
    temporal_knot_spacing : float, optional
        Target spacing between temporal break points for the base model, in
        minutes. The actual spacing is shortened so both ends of the image
        sequence are break points. Default is 120 minutes.
    directional_temporal_knot_spacing : float, optional
        Target temporal break-point spacing for the directional component, in
        minutes. If omitted, it uses `temporal_knot_spacing`. This option is
        only available for directional spatial models.
    minimum_frame_quality : {0, 1, 2, None}, optional
        Minimum `frame_quality` retained in the fit. The default is 2
        (science-ready frames). One also includes detector-valid frames, zero
        includes all classified frames, and None disables filtering. Datasets
        without `frame_quality` are left unchanged.
    tOrder : int, optional
        Order of the temporal spline fit. The default is 2.
    spatial_model : {'legacy', 'legacy_x_azimuth', 'legacy_sh'}, optional
        `legacy` fits the ordinary one-dimensional background. The
        `legacy_x_azimuth` model adds a signed solar-view correction in the
        same legacy coordinate. `legacy_sh` instead adds a coestimated,
        non-axisymmetric real spherical-harmonic correction. Default is
        `legacy`.
    directional_x_knots : array like, optional
        Spatial break points for the signed directional spline. By default it
        uses the ordinary unpadded `sKnots`.
    sh_degree : int, optional
        Maximum spherical-harmonic degree for `legacy_sh`. Default is 4.
    sh_order : int, optional
        Maximum spherical-harmonic order for `legacy_sh`. Default is 4. Only
        terms with m >= 1 and even n-m are retained.
    spencer_taper : bool, optional
        If True, multiply the directional term by the normalized Spencer
        spherical-Chapman profile. The base model is unchanged. Default is
        False.
    variance_weighted : bool, optional
        If True, multiply each fitted row by
        `1/sqrt(imgs['img_variance'])`. Default is False.
    pointing_correction : bool, optional
        Apply shared metadata-only pointing correction before fitting. Default
        False preserves existing behavior. True requires a mapping containing
        WIC, SI12 and SI13 Datasets loaded with pointing_metadata=True.
        Returns a sensor-keyed mapping of fits. Single-sensor correction is
        rejected because shared spacecraft repairs must reach all sensors.
    pointing_options : dict, optional
        Options passed to correct_pointing when pointing_correction=True.
    sigma_model_coefficients : array like, optional
        Supplied nonnegative [a, b] for the empirical residual-scale relation
        `sigma**2 = a*dgmodel**2 + b`. Reevaluate this relation after each fit
        instead of estimating it from bins. Robust weights still start fresh.
        This scale is separate from detector variance and fit uncertainty.
    store_sigma_model : bool, optional
        Save the estimated [a, b] in dgmodel's `sigma_model_coefficients`
        attribute. Default is False. Supplied coefficients are always saved.
    dgsigma_freeze_iteration : {1, 2}, optional
        Default1 estimates sigma after the first fit. Value2 estimates it
        again after the second fit, using the incoming Tukey weights, then
        freezes that sigma field. At least one subsequent fit uses the new
        scale. Incompatible with supplied sigma_model_coefficients.
    initial_robust_weights : array like or xarray.DataArray, optional
        Initial Tukey objective weights in [0, 1], with the same shape as img.
        For a single frame with supplied sigma. Used for the first solve only;
        every fitted pixel, including initial zeros, is reweighted afterwards.
        Numerical damping retains its all-pixels initial scale so only the
        robust initialization changes. Cannot combine with background_prior.
    background_prior : array like or xarray.DataArray, optional
        Preferred background prediction with the same shape as img. Available
        for a single frame with supplied sigma_model_coefficients. Its penalty
        uses fixed bin-population and measurement weights, without Tukey
        rejection. This is empirical shape regularization, not new data.
    background_prior_strength : float, optional
        Nonnegative relative strength of the prediction penalty. Zero disables
        it. One matches the total un-rejected data quadratic. The standard
        numerical damping is centered on the projected prior when enabled.
    store_spatial_coefficients : bool, optional
        Save per-frame spatial spline amplitudes in dgmodel_spatial_coefficients
        for evaluating the same fitted function at another geometry. Default
        is False; this does not change the fitted numerical fields.
    monotonic : bool, optional
        If True, constrain the fitted background to be non-decreasing with
        the legacy coordinate `cos(SZA)/cos(DZA)`. For the directional model,
        this constraint applies only to the base component. Default is False.
    monotonic_fallback : bool, optional
        If True, a failed monotonic solve uses the ordinary solver for that
        IRLS iteration and all remaining iterations. If False, the function
        returns an all-NaN background product with a failure flag. Neither
        case terminates a long processing run. Default is False.
    Returns
    -------
    imgs : xarray.Dataset
        A copy of the input Dataset with new data variables:
            - imgs['dgmodel'] is the B-spline based dayglow model
            - imgs['dgimg'] is the dayglow-corrected image (dayglow subtracked from the input image)
            - imgs['dgweight'] is the weights after the final iteration
            - imgs['dgsigma'] is the modelled spread
            - imgs['dgnorm_residual'] is the residual norm
            - imgs['dgnorm_model'] is the model norm
            - imgs['dgsolver_flag'] is the scalar spatial-solver outcome code
              (0 ordinary, 1 monotonic, 2 ordinary fallback, 3 failed without
              fallback, or 4 failed ordinary fallback)
    '''

    # Shared correction precedes every sensor fit, without changing defaults.
    pointing_correction = kwargs.pop('pointing_correction', False)
    pointing_options = kwargs.pop('pointing_options', None)
    if not isinstance(pointing_correction, (bool, np.bool_)):
        raise ValueError('pointing_correction must be boolean')
    if pointing_options is not None and not isinstance(pointing_options, dict):
        raise ValueError('pointing_options must be a dictionary')
    if pointing_options and not pointing_correction:
        raise ValueError('pointing_options requires pointing_correction=True')
    if isinstance(imgs, Mapping) and not isinstance(imgs, xr.Dataset):
        if not imgs:
            raise ValueError('At least one sensor dataset is required')
        sensors = imgs
        if pointing_correction:
            from fuvpy.src.pointing import correct_pointing
            sensors = correct_pointing(sensors, **(pointing_options or {}))
        return {sensor: backgroundmodel_BS(data, **kwargs) for sensor, data in sensors.items()}
    if pointing_correction:
        raise ValueError('Shared pointing correction requires a WIC/SI12/SI13 mapping')

    # Set keyword arguments to input or default values
    sh_options_given = 'sh_degree' in kwargs or 'sh_order' in kwargs
    sOrder = kwargs.pop('sOrder') if 'sOrder' in kwargs.keys() else 3
    dampingVal = kwargs.pop('dampingVal') if 'dampingVal' in kwargs.keys() else 0
    tukeyVal = kwargs.pop('tukeyVal') if 'tukeyVal' in kwargs.keys() else 5
    stop = kwargs.pop('stop') if 'stop' in kwargs.keys() else 1e-3
    sKnots = kwargs.pop('sKnots') if 'sKnots' in kwargs.keys() else None
    temporal_knot_spacing = kwargs.pop('temporal_knot_spacing') if 'temporal_knot_spacing' in kwargs.keys() else 120
    directional_temporal_knot_spacing = kwargs.pop('directional_temporal_knot_spacing') if 'directional_temporal_knot_spacing' in kwargs.keys() else None
    minimum_frame_quality = kwargs.pop('minimum_frame_quality') if 'minimum_frame_quality' in kwargs.keys() else 2
    tOrder = kwargs.pop('tOrder') if 'tOrder' in kwargs.keys() else 2
    spatial_model = kwargs.pop('spatial_model') if 'spatial_model' in kwargs.keys() else 'legacy'
    directional_x_knots = kwargs.pop('directional_x_knots') if 'directional_x_knots' in kwargs.keys() else None
    sh_degree = kwargs.pop('sh_degree') if 'sh_degree' in kwargs.keys() else 4
    sh_order = kwargs.pop('sh_order') if 'sh_order' in kwargs.keys() else 4
    spencer_taper = bool(kwargs.pop('spencer_taper')) if 'spencer_taper' in kwargs.keys() else False
    variance_weighted = bool(kwargs.pop('variance_weighted')) if 'variance_weighted' in kwargs.keys() else False
    sigma_model_coefficients = kwargs.pop('sigma_model_coefficients', None)
    store_sigma_model = bool(kwargs.pop('store_sigma_model', False))
    dgsigma_freeze_iteration = kwargs.pop('dgsigma_freeze_iteration', 1)
    initial_robust_weights = kwargs.pop('initial_robust_weights', None)
    background_prior = kwargs.pop('background_prior', None)
    background_prior_strength = kwargs.pop('background_prior_strength', 0)
    store_spatial_coefficients = bool(kwargs.pop('store_spatial_coefficients', False))
    monotonic = bool(kwargs.pop('monotonic')) if 'monotonic' in kwargs.keys() else False
    monotonic_fallback = bool(kwargs.pop('monotonic_fallback')) if 'monotonic_fallback' in kwargs.keys() else False
    if kwargs:
        names = ', '.join(sorted(kwargs))
        raise TypeError(f"Unexpected backgroundmodel_BS keyword argument(s): {names}")

    if isinstance(dgsigma_freeze_iteration, (bool, np.bool_)) or not isinstance(dgsigma_freeze_iteration, (int, np.integer)) or dgsigma_freeze_iteration not in [1, 2]:
        raise ValueError("dgsigma_freeze_iteration must be1 or2")
    if dgsigma_freeze_iteration == 2 and sigma_model_coefficients is not None:
        raise ValueError("dgsigma_freeze_iteration=2 cannot reestimate supplied sigma_model_coefficients")
    if sigma_model_coefficients is not None:
        sigma_model_coefficients = np.asarray(sigma_model_coefficients, dtype=float)
        if sigma_model_coefficients.shape != (2,) or np.any(~np.isfinite(sigma_model_coefficients)) or np.any(sigma_model_coefficients < 0) or not np.any(sigma_model_coefficients > 0):
            raise ValueError("sigma_model_coefficients must contain two finite nonnegative values, with at least one positive")
        sigma_model_coefficients = sigma_model_coefficients.copy()
    stored_sigma_model = sigma_model_coefficients
    if isinstance(background_prior_strength, (bool, np.bool_)) or not np.isscalar(background_prior_strength) or not np.isfinite(background_prior_strength) or background_prior_strength < 0:
        raise ValueError("background_prior_strength must be a finite nonnegative scalar")
    use_prior = background_prior_strength > 0
    if use_prior and (background_prior is None or sigma_model_coefficients is None):
        raise ValueError("A background prior requires background_prior and supplied sigma_model_coefficients")
    if initial_robust_weights is not None and (use_prior or sigma_model_coefficients is None):
        raise ValueError("initial_robust_weights requires supplied sigma and cannot combine with a background prior")

    if spatial_model not in ['legacy', 'legacy_x_azimuth', 'legacy_sh']:
        raise ValueError("spatial_model must be 'legacy', 'legacy_x_azimuth', or 'legacy_sh'")

    x_directional_model = spatial_model == 'legacy_x_azimuth'
    sh_directional_model = spatial_model == 'legacy_sh'
    directional_model = x_directional_model or sh_directional_model
    if directional_x_knots is not None and not x_directional_model:
        raise ValueError("directional_x_knots is only available for spatial_model='legacy_x_azimuth'")
    if sh_options_given and not sh_directional_model:
        raise ValueError("sh_degree and sh_order are only available for spatial_model='legacy_sh'")
    if spencer_taper and not directional_model:
        raise ValueError("spencer_taper=True requires a directional spatial model")
    if directional_temporal_knot_spacing is not None and not directional_model:
        raise ValueError("directional_temporal_knot_spacing requires a directional spatial model")
    if minimum_frame_quality is not None and (isinstance(minimum_frame_quality, (bool, np.bool_)) or not isinstance(minimum_frame_quality, (int, np.integer)) or minimum_frame_quality not in [0, 1, 2]):
        raise ValueError("minimum_frame_quality must be 0, 1, 2, or None")

    imgs = imgs.copy(deep=True)

    # Add temporal dimension if missing
    if 'date' not in imgs.img.dims: imgs = imgs.expand_dims('date')

    # Frame quality is assigned by read_idl before WIC reflattening. Subset
    # before constructing the temporal basis so rejected endpoints do not
    # change the fit duration or knot placement.
    frame_quality_filter = 'disabled' if minimum_frame_quality is None else 'not_available'
    if minimum_frame_quality is not None and 'frame_quality' in imgs:
        if imgs['frame_quality'].dims != ('date',):
            raise ValueError("frame_quality must have only the date dimension")
        keep = imgs['frame_quality'].values >= minimum_frame_quality
        if not np.any(keep):
            raise ValueError(f"No frames satisfy minimum_frame_quality={minimum_frame_quality}")
        imgs = imgs.isel(date=keep)
        frame_quality_filter = f"frame_quality >= {minimum_frame_quality}"

    # Spatial knots and viewing angle correction
    if imgs['id'] in ['WIC','SI12','SI13','UVI']:
        sza = imgs['sza'].stack(z=('row','col')).values
        fraction = np.cos(np.deg2rad(sza))/np.cos(np.deg2rad(imgs['dza'].stack(z=('row','col')).values))
        if sKnots is None:
            sKnots = [-3.5,-0.25,0,0.25,1.5,3.5]
    elif imgs['id'] == 'VIS':
        if directional_model:
            raise ValueError("Directional spatial models require IMAGE/FUV viewing geometry")
        sza = imgs['sza'].stack(z=('row','col')).values
        fraction = np.cos(np.deg2rad(sza))
        if sKnots is None: sKnots= [-1,-0.1,0,0.1,0.4,1]
        if monotonic:
            raise ValueError("monotonic fitting requires the legacy FUV coordinate")
    base_breaks = np.asarray(sKnots, dtype=float)
    sKnots = np.r_[np.repeat(base_breaks[0],sOrder),base_breaks, np.repeat(base_breaks[-1],sOrder)]

    sh_keys = []
    if x_directional_model:
        if 'solar_view_sine_factor' not in imgs:
            raise ValueError("spatial_model='legacy_x_azimuth' requires solar_view_sine_factor; load IDL files with viewing_geometry=True")

        directional_x_breaks = base_breaks.copy() if directional_x_knots is None else np.asarray(directional_x_knots, dtype=float)
        if directional_x_breaks.ndim != 1 or directional_x_breaks.size < 2 or np.any(~np.isfinite(directional_x_breaks)) or np.any(np.diff(directional_x_breaks) <= 0):
            raise ValueError("directional_x_knots must be a finite, strictly increasing one-dimensional array")
        directional_x_knots = np.r_[np.repeat(directional_x_breaks[0],sOrder),directional_x_breaks, np.repeat(directional_x_breaks[-1],sOrder)]
        solar_view_sine = imgs['solar_view_sine_factor'].stack(z=('row','col')).values
        directional_taper = _spencer_taper_weights(sza) if spencer_taper else np.ones_like(sza)
    elif sh_directional_model:
        if 'glat' not in imgs or 'glon' not in imgs:
            raise ValueError("spatial_model='legacy_sh' requires glat and glon")

        sh_keys = _nonaxisymmetric_sh_keys(sh_degree, sh_order)
        geographic_latitude = imgs['glat'].stack(z=('row','col')).values
        geographic_longitude = imgs['glon'].stack(z=('row','col')).values
        _, subsolar_longitude = subsol(imgs['date'].values)
        subsolar_longitude = np.asarray(subsolar_longitude)
        solar_relative_longitude = (geographic_longitude - subsolar_longitude[:, None] + 180) % 360 - 180
        directional_taper = _spencer_taper_weights(sza) if spencer_taper else np.ones_like(sza)

    # Minuets since first image
    time=(imgs['date'].values-imgs['date'].values[0])/ np.timedelta64(1, 'm')

    # temporal and spatial size
    n_t = fraction.shape[0]
    n_s = fraction.shape[1]

    # Temporal bases. The directional component may vary on a different time
    # scale while retaining the same spatial support and IRLS weights.
    temporal_breaks, M_base, n_base_tcp, base_tOrder, actual_temporal_knot_spacing = _temporal_spline_basis(time, temporal_knot_spacing, tOrder)
    if directional_model:
        if directional_temporal_knot_spacing is None:
            directional_temporal_knot_spacing = temporal_knot_spacing
        directional_temporal_breaks, M_directional, n_directional_tcp, directional_tOrder, actual_directional_temporal_knot_spacing = _temporal_spline_basis(time, directional_temporal_knot_spacing, tOrder)
    else:
        directional_temporal_breaks = None
        M_directional = None
        n_directional_tcp = 0
        directional_tOrder = None
        actual_directional_temporal_knot_spacing = None

    # Number of spatial control points
    n_base_scp = len(sKnots)-sOrder-1
    if x_directional_model:
        n_directional_scp = len(directional_x_knots)-sOrder-1
    elif sh_directional_model:
        n_directional_scp = 2*len(sh_keys)
    else:
        n_directional_scp = 0
    n_base_coefficients = n_base_scp*n_base_tcp
    n_directional_coefficients = n_directional_scp*n_directional_tcp
    n_coefficients = n_base_coefficients + n_directional_coefficients

    # Index of data to use in model
    fit_support = (sza >= 0) & (np.isfinite(imgs['img'].stack(z=('row','col')).values)) & (fraction>sKnots[0]) & (fraction<sKnots[-1])
    if x_directional_model:
        fit_support &= np.isfinite(solar_view_sine) & np.isfinite(directional_taper)
        fit_support &= (fraction>directional_x_knots[0]) & (fraction<directional_x_knots[-1])
    elif sh_directional_model:
        fit_support &= np.isfinite(geographic_latitude) & np.isfinite(solar_relative_longitude) & np.isfinite(directional_taper)

    monotonic_transform = None
    monotonic_lower = None
    if monotonic and directional_model:
        monotonic_transform, monotonic_lower = _monotonic_base_transform(n_base_scp, n_base_tcp, n_directional_coefficients)

    # Give every populated geometry bin equal total squared influence.
    spatial_bins = np.arange(sKnots[0], sKnots[-1] + 0.1, 0.1)
    ws = _bin_population_weights(fraction, fit_support, spatial_bins)

    # Make everything flat
    d_s = imgs['img'].stack(z=('row','col')).values.flatten()
    ind = fit_support.flatten()
    ws = ws.flatten()
    w = np.ones_like(d_s)

    # Measurement precision is a separate row multiplier. It affects only the
    # background inversion, not the empirical dgsigma regression or the Tukey
    # residual scaling below.
    measurement_w = np.ones_like(d_s)
    if variance_weighted:
        if 'img_variance' not in imgs:
            raise ValueError("variance_weighted=True requires imgs['img_variance']")

        try:
            measurement_variance = imgs['img_variance'].transpose('date', 'row', 'col').stack(z=('row', 'col')).values.flatten()
        except ValueError:
            raise ValueError("img_variance must have date, row, and col dimensions")

        if measurement_variance.shape != d_s.shape:
            raise ValueError("img_variance must have the same shape as img")
        if np.any(~np.isfinite(measurement_variance[ind])):
            raise ValueError("img_variance must be finite at fitted pixels")
        if np.any(measurement_variance[ind] <= 0):
            raise ValueError("img_variance must be strictly positive at fitted pixels")

        measurement_w[ind] = 1/np.sqrt(measurement_variance[ind])

    # Evaluate spatial bases only where image rows contribute to the fit.
    # The temporal basis is applied while accumulating the normal system, so
    # no pixels-by-coefficients whole-orbit matrix is materialized.
    frame_support = [np.flatnonzero(fit_support[frame]) for frame in range(n_t)]
    spatial_designs = []
    for frame, pixel_index in enumerate(frame_support):
        G_base = BSpline(sKnots, np.eye(n_base_scp), sOrder)(fraction[frame, pixel_index])
        if x_directional_model:
            G_directional = BSpline(directional_x_knots, np.eye(n_directional_scp), sOrder)(fraction[frame, pixel_index])
            G_sine = G_directional*(directional_taper[frame, pixel_index]*solar_view_sine[frame, pixel_index])[:, None]
            spatial_designs.append((G_base, G_sine))
        elif sh_directional_model:
            G_sh, _ = _real_sh_design(geographic_latitude[frame, pixel_index], solar_relative_longitude[frame, pixel_index], sh_degree, sh_order)
            G_sh *= directional_taper[frame, pixel_index, None]
            spatial_designs.append((G_base, G_sh))
        else:
            spatial_designs.append((G_base, None))

    if initial_robust_weights is not None:
        if n_t != 1:
            raise ValueError("initial_robust_weights is supported only for a single frame")
        if isinstance(initial_robust_weights, xr.DataArray):
            if initial_robust_weights.dims != imgs.img.dims:
                raise ValueError("initial_robust_weights dimensions must match img")
            xr.align(imgs.img, initial_robust_weights, join='exact')
        initial = np.asarray(initial_robust_weights, dtype=float)
        if initial.shape != imgs.img.shape:
            raise ValueError("initial_robust_weights must match img shape")
        if np.any(~np.isfinite(initial.reshape(-1)[ind])) or np.any((initial.reshape(-1)[ind] < 0) | (initial.reshape(-1)[ind] > 1)) or not np.any(initial.reshape(-1)[ind] > 0):
            raise ValueError("initial_robust_weights must be finite in [0, 1] at fitted pixels, with positive support")
        unweighted_matrix, _ = _accumulate_normal_system(spatial_designs, frame_support, M_base, M_directional, d_s, w, ws, measurement_w, n_s)
        initial_damping_scale = _normal_matrix_scale(unweighted_matrix)
        w[ind] = initial.reshape(-1)[ind]

    # Anchor predictions over this frame's fitted support. Fix the prior
    # quadratic before robust rejection; aurora cannot erase this constraint.
    if use_prior:
        if n_t != 1:
            raise ValueError("background_prior is supported only for a single frame")
        if isinstance(background_prior, xr.DataArray):
            if background_prior.dims != imgs.img.dims:
                raise ValueError("background_prior dimensions must match img")
            xr.align(imgs.img, background_prior, join='exact')
        prior = np.asarray(background_prior, dtype=float)
        if prior.shape != imgs.img.shape or np.any(~np.isfinite(prior.reshape(-1)[ind])):
            raise ValueError("background_prior must match img shape and be finite at fitted pixels")
        prior_matrix, prior_data = _accumulate_normal_system(spatial_designs, frame_support, M_base, M_directional, prior.reshape(-1), w, ws, measurement_w, n_s)
        prior_coefficients = lstsq(prior_matrix, prior_data, lapack_driver='gelsy')[0]
        a, b = sigma_model_coefficients
        prior_sigma = np.sqrt(a*prior.reshape(-1)[ind]**2 + b)
        if np.any(~np.isfinite(prior_sigma)) or np.any(prior_sigma <= 0):
            raise ValueError("Supplied sigma model must give positive finite scales at the prior")
        prior_residual = d_s[ind] - prior.reshape(-1)[ind]
        w[ind] = (1 - np.minimum(1, (prior_residual/(tukeyVal*prior_sigma))**2))**2

    # Iterativaly solve the inverse problem
    diff = 1e10
    iteration = 0
    m = None
    R = None
    use_monotonic = monotonic
    solver_flag = 1 if monotonic else 0
    solver_message = 'success'
    fit_failed = False
    while (diff>stop)&(iteration < 100):
        # Tukey values are objective weights. Spatial and measurement terms
        # are row multipliers. Accumulate their normal system frame by frame
        # so memory follows the input images rather than a much larger
        # pixels-by-coefficients design matrix.
        GTG, GTd = _accumulate_normal_system(spatial_designs, frame_support, M_base, M_directional, d_s, w, ws, measurement_w, n_s)

        # Fix the damping scale before later IRLS iterations change w.
        if R is None:
            gtg_magnitude = initial_damping_scale if initial_robust_weights is not None else _normal_matrix_scale(prior_matrix if use_prior else GTG)
            R = dampingVal*gtg_magnitude*np.eye(n_coefficients)

        if use_prior:
            prior_trace_fraction = background_prior_strength*np.trace(prior_matrix)/(np.trace(GTG) + background_prior_strength*np.trace(prior_matrix))
            GTG += background_prior_strength*prior_matrix
            GTd += background_prior_strength*prior_data + R@prior_coefficients

        if use_monotonic:
            try:
                if use_prior and iteration == 0:
                    # Choose a feasible center for unobserved coefficients.
                    # Keep this solve inside the normal fallback handling.
                    projected_prior = _solve_monotonic_bs(prior_matrix, prior_data, np.zeros_like(prior_matrix), n_base_scp, n_base_tcp, 'normal_system', monotonic_transform, monotonic_lower)
                    GTd += R@(projected_prior - prior_coefficients)
                    prior_coefficients = projected_prior
                m_s = _solve_monotonic_bs(GTG, GTd, R, n_base_scp, n_base_tcp, 'normal_system', monotonic_transform, monotonic_lower)
            except RuntimeError as error:
                if not monotonic_fallback:
                    warnings.warn(f"Monotonic B-spline solve failed; returning an all-NaN background product: {error}", RuntimeWarning, stacklevel=2)
                    solver_flag = 3
                    solver_message = str(error)
                    fit_failed = True
                    break
                warnings.warn(f"Monotonic B-spline solve failed; using the ordinary solver for this and all remaining IRLS iterations: {error}", RuntimeWarning, stacklevel=2)
                use_monotonic = False
                solver_flag = 2
                solver_message = str(error)
        if not use_monotonic:
            try:
                m_s = lstsq(GTG + R, GTd, lapack_driver='gelsy')[0]
            except (np.linalg.LinAlgError, ValueError) as error:
                if solver_flag != 2:
                    raise
                warnings.warn(f"Ordinary fallback solve also failed; returning an all-NaN background product: {error}", RuntimeWarning, stacklevel=2)
                solver_flag = 4
                solver_message = str(error)
                fit_failed = True
                break
        base_coefficients = m_s[:n_base_coefficients].reshape((n_base_scp, n_base_tcp))
        base_amplitude = M_base@base_coefficients.T
        if directional_model:
            directional_coefficients = m_s[n_base_coefficients:].reshape((n_directional_scp, n_directional_tcp))
            directional_amplitude = M_directional@directional_coefficients.T
            mNew = np.hstack([base_amplitude, directional_amplitude])
        else:
            directional_amplitude = None
            mNew = base_amplitude
        dm = np.full((n_t, n_s), np.nan)
        if directional_model:
            dm_base = np.full((n_t, n_s), np.nan)
            dm_directional = np.full((n_t, n_s), np.nan)

        for frame, ((G_base, G_directional), pixel_index) in enumerate(zip(spatial_designs, frame_support)):
            frame_base = G_base@base_amplitude[frame]
            if directional_model:
                frame_directional = G_directional@directional_amplitude[frame]
                dm[frame, pixel_index] = frame_base + frame_directional
                dm_base[frame, pixel_index] = frame_base
                dm_directional[frame, pixel_index] = frame_directional
            else:
                dm[frame, pixel_index] = frame_base

        residuals = (d_s[ind] - dm.flatten()[ind])

        if iteration == 0:
            sigma = np.full_like(d_s,np.nan)
        if sigma_model_coefficients is None and iteration < dgsigma_freeze_iteration:
            if store_sigma_model:
                sigma[ind], stored_sigma_model = _noisemodel_binned(fraction.flatten()[ind], d_s[ind], dm.flatten()[ind], w[ind]*ws[ind]**2, sKnots, return_coefficients=True)
            else:
                sigma[ind] = _noisemodel_binned(fraction.flatten()[ind], d_s[ind], dm.flatten()[ind], w[ind]*ws[ind]**2, sKnots)
        if sigma_model_coefficients is not None:
            a, b = sigma_model_coefficients
            sigma[ind] = np.sqrt(a*dm.flatten()[ind]**2 + b)
            if np.any(~np.isfinite(sigma[ind])) or np.any(sigma[ind] <= 0):
                raise ValueError("Supplied sigma model must give positive finite scales at fitted pixels")
        w[ind] = (1 - np.minimum(1,(residuals/(tukeyVal*sigma[ind]))**2))**2

        if m is not None: diff = np.sqrt(np.mean((mNew-m)**2))/(1+np.sqrt(np.mean(mNew**2)))
        if sigma_model_coefficients is None and iteration < dgsigma_freeze_iteration:
            # Consume the final sigma estimate in another weighted solve.
            diff = 1e10
        m = mNew
        iteration += 1

    if iteration == 100: print('Warning: Model did not converge. Truncated after 100 iterations')

    # A failed constrained fit must not leave a plausible-looking product.
    # NaNs preserve the output schema while dgsolver_flag records why the fit
    # is unavailable. The long-running caller can save the failure and move on.
    if fit_failed:
        shape = (n_t, len(imgs.row), len(imgs.col))
        failed_field = np.full(shape, np.nan)
        imgs['dgmodel'] = (['date','row','col'], failed_field.copy())
        imgs['dgimg'] = (['date','row','col'], failed_field.copy())
        imgs['dgweight'] = (['date','row','col'], failed_field.copy())
        imgs['dgsigma'] = (['date','row','col'], failed_field.copy())
        imgs['dgnorm_residual'] = np.nan
        imgs['dgnorm_model'] = np.nan
        if directional_model:
            imgs['dgmodel_base'] = (['date','row','col'], failed_field.copy())
            if x_directional_model:
                imgs['dgmodel_sine'] = (['date','row','col'], failed_field.copy())
            else:
                imgs['dgmodel_sh'] = (['date','row','col'], failed_field.copy())
    else:
        normM = np.sqrt(np.average(m_s**2))
        normR = np.nan if use_prior and not np.any(w[ind] > 0) else np.sqrt(np.average(residuals**2,weights=w[ind]))
        imgs['dgmodel'] = (['date','row','col'],(dm).reshape((n_t,len(imgs.row),len(imgs.col))))
        imgs['dgimg'] = imgs['img']-imgs['dgmodel']
        imgs['dgweight'] = (['date','row','col'],(w).reshape((n_t,len(imgs.row),len(imgs.col))))
        imgs['dgsigma'] = (['date','row','col'],(sigma).reshape((n_t,len(imgs.row),len(imgs.col))))
        imgs['dgnorm_residual'] = normR
        imgs['dgnorm_model'] = normM
        if directional_model:
            imgs['dgmodel_base'] = (['date','row','col'],dm_base.reshape((n_t,len(imgs.row),len(imgs.col))))
            if x_directional_model:
                imgs['dgmodel_sine'] = (['date','row','col'],dm_directional.reshape((n_t,len(imgs.row),len(imgs.col))))
            else:
                imgs['dgmodel_sh'] = (['date','row','col'],dm_directional.reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['dgsolver_flag'] = solver_flag
    if directional_model and spencer_taper:
        imgs['spencer_taper_factor'] = (['date','row','col'],directional_taper.reshape((n_t,len(imgs.row),len(imgs.col))))
    fit_converged = False if fit_failed else iteration < 100 and np.isfinite(diff) and diff <= stop
    imgs['dgfit_iterations'] = iteration
    imgs['dgfit_converged'] = fit_converged
    if store_spatial_coefficients:
        amplitudes = np.full((n_t, n_base_scp + n_directional_scp), np.nan) if fit_failed else mNew
        imgs['dgmodel_spatial_coefficients'] = (['date', 'dg_spatial_coefficient'], amplitudes)
        imgs['dgmodel_spatial_coefficients'].attrs = {'long_name': 'Per-frame spatial amplitudes; base first, then directional', 'base_count': n_base_scp, 'directional_count': n_directional_scp}

    # Do not return extrapolated models or fitting weights outside the pixels
    # that actually contributed to the inverse problem.
    fit_support = xr.DataArray(fit_support.reshape((n_t, len(imgs.row), len(imgs.col))), dims=['date', 'row', 'col'], coords={'date': imgs.date, 'row': imgs.row, 'col': imgs.col})
    imgs['dgmodel'] = xr.where(fit_support,imgs['dgmodel'],np.nan)
    imgs['dgimg'] = xr.where(fit_support,imgs['dgimg'],np.nan)
    imgs['dgweight'] = xr.where(fit_support,imgs['dgweight'],np.nan)
    imgs['dgsigma'] = xr.where(fit_support,imgs['dgsigma'],np.nan)
    if directional_model:
        imgs['dgmodel_base'] = xr.where(fit_support,imgs['dgmodel_base'],np.nan)
        if x_directional_model:
            imgs['dgmodel_sine'] = xr.where(fit_support,imgs['dgmodel_sine'],np.nan)
        else:
            imgs['dgmodel_sh'] = xr.where(fit_support,imgs['dgmodel_sh'],np.nan)

    # Add attributes
    monotonic_constraint = 'none'
    if monotonic and solver_flag == 1:
        monotonic_constraint = 'legacy_base_nondecreasing' if directional_model else 'legacy_coordinate_nondecreasing'
    elif solver_flag in [3, 4]:
        monotonic_constraint = 'not_available'
    imgs['dgmodel'].attrs = {
        'long_name': 'BS model',
        'spatial_model': spatial_model,
        'spatial_knots': base_breaks,
        'spatial_order': int(sOrder),
        'temporal_knots_minutes_from_first_image': temporal_breaks,
        'temporal_knot_count': int(len(temporal_breaks)),
        'temporal_order': int(base_tOrder),
        'temporal_knot_spacing_target_minutes': float(temporal_knot_spacing),
        'temporal_knot_spacing_actual_minutes': float(actual_temporal_knot_spacing),
        'frame_quality_filter': frame_quality_filter,
        'noise_model': 'binned' if sigma_model_coefficients is None else 'supplied_quadratic',
        'spatial_weighting': 'bin_population',
        'dgsigma_estimation': 'initial_fit' if sigma_model_coefficients is None else 'supplied_relation_each_iteration',
        'robust_weighting': 'tukey_bisquare',
        'tukey_value': float(tukeyVal),
        'measurement_variance': 'img_variance' if variance_weighted else 'not_used',
        'variance_weighted': 'true' if variance_weighted else 'false',
        'monotonic_requested': 'true' if monotonic else 'false',
        'monotonic_constraint': monotonic_constraint,
        'monotonic_solver': 'normal' if monotonic else 'not_used',
        'monotonic_fallback': 'true' if monotonic_fallback else 'false',
        'spencer_taper': 'enabled' if spencer_taper else 'disabled',
        'damping_value': float(dampingVal),
        'damping_scale': float(gtg_magnitude),
        'convergence_tolerance': float(stop),
        'maximum_iterations': 100,
        'fit_iterations': int(iteration),
        'fit_converged': int(fit_converged)
    }
    if stored_sigma_model is not None:
        imgs['dgmodel'].attrs['sigma_model_coefficients'] = stored_sigma_model
        imgs['dgmodel'].attrs['sigma_model_equation'] = 'residual_scale_squared = a*background_squared + b; not detector variance or background uncertainty'
    if dgsigma_freeze_iteration == 2:
        imgs['dgmodel'].attrs.update({'dgsigma_estimation': 'first_two_fits_then_frozen', 'dgsigma_freeze_iteration': 2})
    if initial_robust_weights is not None:
        imgs['dgmodel'].attrs['robust_initialization'] = 'supplied initial objective weights; all fitted pixels reweighted after every solve'
    if use_prior:
        imgs['dgmodel'].attrs.update({'background_prior_strength': float(background_prior_strength), 'background_prior_penalty': 'Fixed population/measurement-weighted squared prediction departure; numerical damping centered on projected prior', 'robust_initialization': 'fresh residuals against supplied prior'})
        imgs['dgmodel'].attrs['background_prior_trace_fraction'] = float(prior_trace_fraction)
    if directional_model:
        imgs['dgmodel'].attrs['directional_temporal_knots_minutes_from_first_image'] = directional_temporal_breaks
        imgs['dgmodel'].attrs['directional_temporal_knot_count'] = int(len(directional_temporal_breaks))
        imgs['dgmodel'].attrs['directional_temporal_order'] = int(directional_tOrder)
        imgs['dgmodel'].attrs['directional_temporal_knot_spacing_target_minutes'] = float(directional_temporal_knot_spacing)
        imgs['dgmodel'].attrs['directional_temporal_knot_spacing_actual_minutes'] = float(actual_directional_temporal_knot_spacing)
    imgs['dgimg'].attrs = {'long_name': 'BS corrected image'}
    imgs['dgweight'].attrs = {'long_name': 'BS model weights'}
    imgs['dgsigma'].attrs = {'long_name': 'BS model spread'}
    imgs['dgnorm_residual'].attrs = {'long_name': 'BS model residual norm'}
    imgs['dgnorm_model'].attrs = {'long_name': 'BS model solution norm'}
    imgs['dgsolver_flag'].attrs = {
        'long_name': 'BS spatial solver outcome',
        'flag_values': np.array([0, 1, 2, 3, 4], dtype=np.int8),
        'flag_meanings': 'standard_requested monotonic_succeeded monotonic_failed_standard_fallback monotonic_failed_no_fallback monotonic_failed_standard_fallback_failed',
        'solver_message': solver_message
    }
    if x_directional_model:
        imgs['dgmodel'].attrs['directional_x_knots'] = directional_x_breaks
        imgs['dgmodel_base'].attrs = {'long_name': 'Azimuth-independent legacy-coordinate background component'}
        imgs['dgmodel_sine'].attrs = {'long_name': 'Sine-like solar-view directional background component'}
    elif sh_directional_model:
        terms = [f"cos({degree},{order})" for degree, order in sh_keys] + [f"sin({degree},{order})" for degree, order in sh_keys]
        imgs['dgmodel'].attrs['sh_degree'] = int(sh_degree)
        imgs['dgmodel'].attrs['sh_order'] = int(sh_order)
        imgs['dgmodel'].attrs['sh_terms'] = ' '.join(terms)
        imgs['dgmodel'].attrs['sh_axisymmetric_terms'] = 'excluded'
        imgs['dgmodel'].attrs['sh_coordinate_system'] = 'geographic_latitude_solar_relative_longitude'
        imgs['dgmodel_base'].attrs = {'long_name': 'Azimuth-independent legacy-coordinate background component'}
        imgs['dgmodel_sh'].attrs = {'long_name': 'Spencer-tapered non-axisymmetric spherical-harmonic background component'}
    if directional_model and spencer_taper:
        imgs['spencer_taper_factor'].attrs = {'long_name': 'Normalized spherical-Chapman production factor', 'units': '1'}
    imgs['dgfit_iterations'].attrs = {'long_name': 'Number of BS IRLS iterations'}
    imgs['dgfit_converged'].attrs = {'long_name': 'BS IRLS convergence status'}
    return imgs

def _noise_bins(sKnots):
    """Return one nightside and fourteen 0.25-wide dayside bins."""

    return np.r_[sKnots[0], np.arange(0, 3.5 + 0.25, 0.25)]


def _bin_population_weights(fraction, valid, bins):
    """Equalize populated-bin influence in each image frame.

    These are row multipliers, so the least-squares objective squares them.
    The common per-frame normalization makes their mean square one over fitted
    pixels and therefore preserves the overall scale relative to damping.
    """

    weights = np.full_like(fraction, np.nan)

    for frame in range(fraction.shape[0]):
        frame_fraction = fraction[frame, valid[frame]]
        if frame_fraction.size == 0:
            continue

        count, _, bin_index = binned_statistic(frame_fraction, frame_fraction, statistic='count', bins=bins)

        populated_bins = np.count_nonzero(count)
        bin_count = count[bin_index - 1]
        normalization = frame_fraction.size / populated_bins
        weights[frame, valid[frame]] = np.sqrt(normalization / bin_count)

    return weights


def _binned_noise_statistics(fraction, d, dm, w, sKnots):
    """Calculate one model/residual variance pair per populated bin."""

    bins = _noise_bins(sKnots)
    binnumber = np.digitize(fraction, bins)

    model_squared = []
    residual_squared = []

    for bin_index in range(1, len(bins)):
        in_bin = binnumber == bin_index
        if np.sum(w[in_bin]) == 0:
            continue

        model_squared.append(np.average(dm[in_bin]**2, weights=w[in_bin]))
        residual_squared.append(np.average((d[in_bin] - dm[in_bin])**2, weights=w[in_bin]))

    return np.asarray(model_squared), np.asarray(residual_squared)


def _noisemodel_binned(fraction, d, dm, w, sKnots, return_coefficients=False):
    """Estimate residual spread without repeating bin statistics per pixel."""

    model_squared, residual_squared = _binned_noise_statistics(fraction, d, dm, w, sKnots)

    if model_squared.size < 2:
        raise ValueError("At least two populated geometry bins are required")

    def func(x, a, b):
        return a*x+b

    coefficients, _ = curve_fit(func, model_squared, residual_squared, bounds=(0, np.inf))

    sigma = np.sqrt(func(dm**2, *coefficients))
    return (sigma, coefficients) if return_coefficients else sigma
