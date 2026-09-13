#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:28:25 2022

@author: aohma
"""
from pathlib import Path
import warnings

import numpy as np
import xarray as xr

from scipy.interpolate import BSpline
from scipy.optimize import curve_fit, lsq_linear
from scipy.stats import binned_statistic
from scipy.linalg import lstsq


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


def _monotonic_base_transform(n_base_scp, n_total_scp, n_tcp):
    """Constrain only the leading one-dimensional base coefficients.

    Coefficients after the base belong to the signed directional model and
    remain free.
    """

    n_coefficients = n_total_scp*n_tcp
    transform = np.eye(n_coefficients)
    lower = np.full(n_coefficients, -np.inf)

    for temporal in range(n_tcp):
        baseline = temporal
        for spatial in range(n_base_scp):
            row = spatial*n_tcp + temporal
            transform[row, :] = 0
            transform[row, baseline] = 1

            for increment in range(1, spatial + 1):
                column = increment*n_tcp + temporal
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
    else:
        # Operational calculation: reduce all pixels to the small normal
        # system before applying the same coefficient bounds.
        normal_matrix = weighted_design.T@weighted_design + damping
        normal_data = weighted_design.T@weighted_data

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

def backgroundmodel_BS(imgs,**kwargs):
    '''
    Function to model the FUV dayglow and subtract it from the input image

    Parameters
    ----------
    imgs : xarray.Dataset
        Dataset with the FUV images, imported by readFUVimage()
    inImg : str, optional
        Name of the input image to be used in the model. The default is 'img'.
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
    n_tKnots : int, optional
        Number of temporal knots, equally spaced between the endpoints. The default is 2 (only knots at endpoints)
    tOrder : int, optional
        Order of the temporal spline fit. The default is 2.
    spatial_model : {'legacy', 'legacy_x_azimuth'}, optional
        `legacy` fits the ordinary one-dimensional background. The
        `legacy_x_azimuth` model adds a signed solar-view correction in the
        same legacy coordinate. Default is `legacy`.
    directional_x_knots : array like, optional
        Spatial break points for the signed directional spline. By default it
        uses the ordinary unpadded `sKnots`.
    spencer_taper : bool, optional
        If True, multiply the directional term by the normalized Spencer
        spherical-Chapman profile. The base model is unchanged. Default is
        False.
    noise_model : {'legacy', 'binned'}, optional
        Method used to estimate `dgsigma`. `legacy` reproduces the released
        pixel-expanded regression. `binned` fits one residual variance and
        one mean model-squared value per geometry bin. Default is `legacy`.
    spatial_weighting : {'gaussian', 'bin_population'}, optional
        Spatial row multiplier used by the least-squares fit. `gaussian`
        reproduces the released inverse-Gaussian multiplier.
        `bin_population` gives every populated 0.1-wide legacy-coordinate bin
        equal total squared influence, independently in each frame. Default
        is `gaussian`.
    dgsigma_freeze_iteration : int, optional
        Last zero-based IRLS iteration on which `dgsigma` is recalculated.
        Zero retains the estimate from the initial fit. The default is 1,
        reproducing the released calculation through iteration 1.
    input_variance : str, optional
        Name of the Dataset variable containing measurement variance in the
        input image units squared. Required when `variance_weighted=True`.
    variance_weighted : bool, optional
        If True, multiply each fitted row by
        `1/sqrt(imgs[input_variance])`. Default is False. With nonzero damping,
        reassess `dampingVal` because absolute measurement weights change the
        balance between the data and regularization.
    monotonic : bool, optional
        If True, constrain the fitted background to be non-decreasing with
        the legacy coordinate `cos(SZA)/cos(DZA)`. For the directional model,
        this constraint applies only to the base component. Default is False.
    monotonic_solver : {'normal', 'full'}, optional
        `normal` first reduces the weighted image rows to a small normal
        system and is the operational default. `full` applies the same bounds
        to all weighted rows as a reference calculation.
    monotonic_fallback : bool, optional
        If True, a failed monotonic solve uses the ordinary solver for that
        IRLS iteration and all remaining iterations. If False, the function
        returns an all-NaN background product with a failure flag. Neither
        case terminates a long processing run. Default is False.
    diagnostics : bool, optional
        If True, add scalar `dgfit_iterations` and `dgfit_converged` fields.
        Failed solver products include these fields even when diagnostics is
        False. Default is False.
    inplace : bool, optional
        If True, update the Dataset in place. Default is False (return a new Dataset).


    Returns
    -------
    imgs : xarray.Dataset, if inplace is False
        A copy of the input Dataset with four new data_vars:
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

    # Set keyword arguments to input or default values    
    inImg = kwargs.pop('inImg') if 'inImg' in kwargs.keys() else 'img'
    sOrder = kwargs.pop('sOrder') if 'sOrder' in kwargs.keys() else 3
    dampingVal = kwargs.pop('dampingVal') if 'dampingVal' in kwargs.keys() else 0
    tukeyVal = kwargs.pop('tukeyVal') if 'tukeyVal' in kwargs.keys() else 5
    stop = kwargs.pop('stop') if 'stop' in kwargs.keys() else 1e-3
    sKnots = kwargs.pop('sKnots') if 'sKnots' in kwargs.keys() else None
    n_tKnots = kwargs.pop('n_tKnots') if 'n_tKnots' in kwargs.keys() else 2
    tOrder = kwargs.pop('tOrder') if 'tOrder' in kwargs.keys() else 2
    spatial_model = kwargs.pop('spatial_model') if 'spatial_model' in kwargs.keys() else 'legacy'
    directional_x_knots = kwargs.pop('directional_x_knots') if 'directional_x_knots' in kwargs.keys() else None
    spencer_taper = bool(kwargs.pop('spencer_taper')) if 'spencer_taper' in kwargs.keys() else False
    noise_model = kwargs.pop('noise_model') if 'noise_model' in kwargs.keys() else 'legacy'
    spatial_weighting = kwargs.pop('spatial_weighting') if 'spatial_weighting' in kwargs.keys() else 'gaussian'
    dgsigma_freeze_iteration = kwargs.pop('dgsigma_freeze_iteration') if 'dgsigma_freeze_iteration' in kwargs.keys() else 1
    input_variance = kwargs.pop('input_variance') if 'input_variance' in kwargs.keys() else None
    variance_weighted = bool(kwargs.pop('variance_weighted')) if 'variance_weighted' in kwargs.keys() else False
    monotonic = bool(kwargs.pop('monotonic')) if 'monotonic' in kwargs.keys() else False
    monotonic_solver = kwargs.pop('monotonic_solver') if 'monotonic_solver' in kwargs.keys() else 'normal'
    monotonic_fallback = bool(kwargs.pop('monotonic_fallback')) if 'monotonic_fallback' in kwargs.keys() else False
    diagnostics = bool(kwargs.pop('diagnostics')) if 'diagnostics' in kwargs.keys() else False
    inplace = bool(kwargs.pop('inplace')) if 'inplace' in kwargs.keys() else False

    # These limits belonged to the old fit API. Current files apply detector
    # and viewing-angle masks in read_idl(), but the paper scripts still pass
    # the historical arguments and are kept runnable for reference.
    obsolete = [name for name in ['dzalim', 'minlat'] if name in kwargs]
    for name in obsolete:
        kwargs.pop(name)
    if obsolete:
        warnings.warn(f"{', '.join(obsolete)} no longer restrict backgroundmodel_BS; apply geometry limits when loading the images", DeprecationWarning, stacklevel=2)
    if kwargs:
        names = ', '.join(sorted(kwargs))
        raise TypeError(f"Unexpected backgroundmodel_BS keyword argument(s): {names}")

    if noise_model not in ['legacy', 'binned']:
        raise ValueError("noise_model must be 'legacy' or 'binned'")

    if spatial_weighting not in ['gaussian', 'bin_population']:
        raise ValueError("spatial_weighting must be 'gaussian' or 'bin_population'")

    if spatial_model not in ['legacy', 'legacy_x_azimuth']:
        raise ValueError("spatial_model must be 'legacy' or 'legacy_x_azimuth'")

    directional_model = spatial_model == 'legacy_x_azimuth'
    if directional_x_knots is not None and not directional_model:
        raise ValueError("directional_x_knots is only available for spatial_model='legacy_x_azimuth'")
    if spencer_taper and not directional_model:
        raise ValueError("spencer_taper=True is only available for spatial_model='legacy_x_azimuth'")

    if (isinstance(dgsigma_freeze_iteration, (bool, np.bool_)) or not isinstance(dgsigma_freeze_iteration, (int, np.integer)) or dgsigma_freeze_iteration < 0):
        raise ValueError("dgsigma_freeze_iteration must be a non-negative integer")

    if monotonic_solver not in ['normal', 'full']:
        raise ValueError("monotonic_solver must be 'normal' or 'full'")

    if not inplace: imgs = imgs.copy(deep=True)

    # Add temporal dimension if missing
    if len(imgs.sizes)==2: imgs = imgs.expand_dims('date')

    # Spatial knots and viewing angle correction
    if imgs['id'] in ['WIC','SI12','SI13','UVI']:
        sza = imgs['sza'].stack(z=('row','col')).values
        fraction = np.cos(np.deg2rad(sza))/np.cos(np.deg2rad(imgs['dza'].stack(z=('row','col')).values))
        if sKnots is None:
            sKnots = [-3.5,-0.25,0,0.25,1.5,3.5]
    elif imgs['id'] == 'VIS':
        if directional_model:
            raise ValueError("spatial_model='legacy_x_azimuth' requires IMAGE/FUV viewing geometry")
        sza = imgs['sza'].stack(z=('row','col')).values
        fraction = np.cos(np.deg2rad(sza))
        if sKnots is None: sKnots= [-1,-0.1,0,0.1,0.4,1]
        if spatial_weighting == 'bin_population':
            raise ValueError("bin_population spatial weighting requires the legacy FUV coordinate")
        if monotonic:
            raise ValueError("monotonic fitting requires the legacy FUV coordinate")
    base_breaks = np.asarray(sKnots, dtype=float)
    sKnots = np.r_[np.repeat(base_breaks[0],sOrder),base_breaks, np.repeat(base_breaks[-1],sOrder)]

    if directional_model:
        if 'solar_view_sine_factor' not in imgs:
            raise ValueError("spatial_model='legacy_x_azimuth' requires solar_view_sine_factor; load IDL files with viewing_geometry=True")

        directional_x_breaks = base_breaks.copy() if directional_x_knots is None else np.asarray(directional_x_knots, dtype=float)
        if directional_x_breaks.ndim != 1 or directional_x_breaks.size < 2 or np.any(~np.isfinite(directional_x_breaks)) or np.any(np.diff(directional_x_breaks) <= 0):
            raise ValueError("directional_x_knots must be a finite, strictly increasing one-dimensional array")
        directional_x_knots = np.r_[np.repeat(directional_x_breaks[0],sOrder),directional_x_breaks, np.repeat(directional_x_breaks[-1],sOrder)]
        solar_view_sine = imgs['solar_view_sine_factor'].stack(z=('row','col')).values
        directional_taper = _spencer_taper_weights(sza) if spencer_taper else np.ones_like(sza)

    # Minuets since first image
    time=(imgs['date'].values-imgs['date'].values[0])/ np.timedelta64(1, 'm')

    # temporal and spatial size
    n_t = fraction.shape[0]
    n_s = fraction.shape[1]

    # Temporal knots
    if n_t==1:
        tKnots = np.linspace(time[0], time[-1]+1, n_tKnots)
        tOrder = 0
    else:
        tKnots = np.linspace(time[0], time[-1], n_tKnots)
        tKnots = np.r_[np.repeat(tKnots[0],tOrder),tKnots, np.repeat(tKnots[-1],tOrder)]

    # Number of control points
    n_tcp = len(tKnots)-tOrder-1
    n_base_scp = len(sKnots)-sOrder-1
    n_directional_scp = len(directional_x_knots)-sOrder-1 if directional_model else 0
    n_scp = n_base_scp + n_directional_scp

    # Temporal design matix
    M = BSpline(tKnots, np.eye(n_tcp), tOrder)(time)

    # Spatial design matrix
    G_g=[]
    G_base_g=[]
    G_sine_g=[]
    G_s=[]
    for i in range(n_t):
        G_base = BSpline(sKnots, np.eye(n_base_scp), sOrder)(fraction[i,:])
        if directional_model:
            G_directional = BSpline(directional_x_knots, np.eye(n_directional_scp), sOrder)(fraction[i,:])
            G_sine = G_directional*(directional_taper[i]*solar_view_sine[i])[:, None]
            G = np.hstack([G_base, G_sine])
            G_base_g.append(G_base)
            G_sine_g.append(G_sine)
        else:
            G = G_base
        G_t = np.zeros((n_s, n_scp*n_tcp))
        for j in range(n_tcp):
            G_t[:, np.arange(j, n_scp*n_tcp, n_tcp)] = G*M[i, j]

        G_g.append(G)
        G_s.append(G_t)
    G_s=np.vstack(G_s)

    # Index of data to use in model
    fit_support = (sza >= 0) & (np.isfinite(imgs[inImg].stack(z=('row','col')).values)) & (fraction>sKnots[0]) & (fraction<sKnots[-1])
    if directional_model:
        fit_support &= np.isfinite(solar_view_sine) & np.isfinite(directional_taper)
        fit_support &= (fraction>directional_x_knots[0]) & (fraction<directional_x_knots[-1])

    monotonic_transform = None
    monotonic_lower = None
    if monotonic and directional_model:
        monotonic_transform, monotonic_lower = _monotonic_base_transform(n_base_scp, n_scp, n_tcp)

    # Spatial sampling weights
    if spatial_weighting == 'gaussian':
        ws = 1/np.exp(-0.5*(fraction/1.5)**2)
    else:
        spatial_bins = np.arange(sKnots[0], sKnots[-1] + 0.1, 0.1)
        ws = _bin_population_weights(fraction, fit_support, spatial_bins)

    # Make everything flat
    d_s = imgs[inImg].stack(z=('row','col')).values.flatten()
    ind = fit_support.flatten()
    ws = ws.flatten()
    w = np.ones_like(d_s)

    # Measurement precision is a separate row multiplier. It affects only the
    # background inversion, not the empirical dgsigma regression or the Tukey
    # residual scaling below.
    measurement_w = np.ones_like(d_s)
    if variance_weighted:
        if not isinstance(input_variance, str):
            raise ValueError("input_variance must name a Dataset variable when variance_weighted=True")
        if input_variance not in imgs:
            raise ValueError(f"Input variance variable '{input_variance}' is not in imgs")

        try:
            measurement_variance = imgs[input_variance].transpose('date', 'row', 'col').stack(z=('row', 'col')).values.flatten()
        except ValueError:
            raise ValueError("Input variance must have date, row, and col dimensions")

        if measurement_variance.shape != d_s.shape:
            raise ValueError("Input variance must have the same shape as the image")
        if np.any(~np.isfinite(measurement_variance[ind])):
            raise ValueError("Input variance must be finite at fitted pixels")
        if np.any(measurement_variance[ind] <= 0):
            raise ValueError("Input variance must be strictly positive at fitted pixels")

        measurement_w[ind] = 1/np.sqrt(measurement_variance[ind])

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
        # Keep robust, spatial-sampling, and measurement-precision roles
        # visible before forming the weighted least-squares system. Tukey
        # values are objective weights; the other two are row multipliers.
        fit_w = np.sqrt(w[ind]) * ws[ind] * measurement_w[ind]
        weighted_design = G_s[ind, :] * fit_w[:, None]
        weighted_data = d_s[ind] * fit_w
        GTG = weighted_design.T@weighted_design
        GTd = weighted_design.T@weighted_data

        # Fix the damping scale before later IRLS iterations change w.
        if R is None:
            gtg_magnitude = _normal_matrix_scale(GTG)
            R = dampingVal*gtg_magnitude*np.eye(G_s.shape[1])

        if use_monotonic:
            try:
                m_s = _solve_monotonic_bs(weighted_design, weighted_data, R, n_scp, n_tcp, monotonic_solver, monotonic_transform, monotonic_lower)
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
        mNew    = M@m_s.reshape((n_scp, n_tcp)).T
        dm=[]
        for i, tt in enumerate(time):
            dm.append(G_g[i]@mNew[i, :])

        dm=np.array(dm).squeeze()
        if directional_model:
            dm_base = np.array([G_base_g[i]@mNew[i, :n_base_scp] for i in range(n_t)]).squeeze()
            dm_sine = np.array([G_sine_g[i]@mNew[i, n_base_scp:] for i in range(n_t)]).squeeze()
        residuals = (d_s[ind] - dm.flatten()[ind])

        if iteration <= dgsigma_freeze_iteration:
            sigma = np.full_like(d_s,np.nan)
            if noise_model == 'legacy':
                sigma[ind] = _noisemodel(fraction.flatten()[ind], d_s[ind], dm.flatten()[ind], w[ind]*ws[ind]**2, sKnots)
            else:
                sigma[ind] = _noisemodel_binned(fraction.flatten()[ind], d_s[ind], dm.flatten()[ind], w[ind]*ws[ind]**2, sKnots)
        w[ind] = (1 - np.minimum(1,(residuals/(tukeyVal*sigma[ind]))**2))**2

        if m is not None: diff = np.sqrt(np.mean((mNew-m)**2))/(1+np.sqrt(np.mean(mNew**2)))
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
            imgs['dgmodel_sine'] = (['date','row','col'], failed_field.copy())
    else:
        normM = np.sqrt(np.average(m_s**2))
        normR = np.sqrt(np.average(residuals**2,weights=w[ind]))
        imgs['dgmodel'] = (['date','row','col'],(dm).reshape((n_t,len(imgs.row),len(imgs.col))))
        imgs['dgimg'] = imgs[inImg]-imgs['dgmodel']
        imgs['dgweight'] = (['date','row','col'],(w).reshape((n_t,len(imgs.row),len(imgs.col))))
        imgs['dgsigma'] = (['date','row','col'],(sigma).reshape((n_t,len(imgs.row),len(imgs.col))))
        imgs['dgnorm_residual'] = normR
        imgs['dgnorm_model'] = normM
        if directional_model:
            imgs['dgmodel_base'] = (['date','row','col'],dm_base.reshape((n_t,len(imgs.row),len(imgs.col))))
            imgs['dgmodel_sine'] = (['date','row','col'],dm_sine.reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['dgsolver_flag'] = solver_flag
    if directional_model and spencer_taper:
        imgs['spencer_taper_factor'] = (['date','row','col'],directional_taper.reshape((n_t,len(imgs.row),len(imgs.col))))
    if diagnostics or fit_failed:
        imgs['dgfit_iterations'] = iteration
        imgs['dgfit_converged'] = False if fit_failed else iteration < 100 and np.isfinite(diff) and diff <= stop

    # Do not return extrapolated models or fitting weights outside the pixels
    # that actually contributed to the inverse problem.
    fit_support = xr.DataArray(fit_support.reshape((n_t, len(imgs.row), len(imgs.col))), dims=['date', 'row', 'col'], coords={'date': imgs.date, 'row': imgs.row, 'col': imgs.col})
    imgs['dgmodel'] = xr.where(fit_support,imgs['dgmodel'],np.nan)
    imgs['dgimg'] = xr.where(fit_support,imgs['dgimg'],np.nan)
    imgs['dgweight'] = xr.where(fit_support,imgs['dgweight'],np.nan)
    imgs['dgsigma'] = xr.where(fit_support,imgs['dgsigma'],np.nan)
    if directional_model:
        imgs['dgmodel_base'] = xr.where(fit_support,imgs['dgmodel_base'],np.nan)
        imgs['dgmodel_sine'] = xr.where(fit_support,imgs['dgmodel_sine'],np.nan)

    # Add attributes
    monotonic_constraint = 'none'
    if monotonic and solver_flag == 1:
        monotonic_constraint = 'legacy_base_nondecreasing' if directional_model else 'legacy_coordinate_nondecreasing'
    elif solver_flag in [3, 4]:
        monotonic_constraint = 'not_available'
    imgs['dgmodel'].attrs = {'long_name': 'BS model', 'spatial_model': spatial_model, 'monotonic_constraint': monotonic_constraint}
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
    if directional_model:
        imgs['dgmodel'].attrs['directional_x_knots'] = directional_x_breaks
        imgs['dgmodel'].attrs['spencer_taper'] = 'enabled' if spencer_taper else 'disabled'
        imgs['dgmodel_base'].attrs = {'long_name': 'Azimuth-independent legacy-coordinate background component'}
        imgs['dgmodel_sine'].attrs = {'long_name': 'Sine-like solar-view directional background component'}
        if spencer_taper:
            imgs['spencer_taper_factor'].attrs = {'long_name': 'Normalized spherical-Chapman production factor', 'units': '1'}
    if diagnostics or fit_failed:
        imgs['dgfit_iterations'].attrs = {'long_name': 'Number of BS IRLS iterations'}
        imgs['dgfit_converged'].attrs = {'long_name': 'BS IRLS convergence status'}
    # Return the new DataSet if not inplace 
    if not inplace:
        return imgs

def _noisemodel(fraction,d,dm,w,sKnots):

    # Binned RMSE
    bins = np.r_[sKnots[0],np.arange(0,sKnots[-1]+0.25,0.25)]
    binnumber = np.digitize(fraction,bins)

    rmse=np.full_like(d,np.nan)
    for i in range(1,len(bins)):
        ind = binnumber==i
        if np.sum(w[ind])==0:
            rmse[ind]=np.nan
        else:
            rmse[ind]=np.sqrt(np.average((d[ind]-dm[ind])**2,weights=w[ind]))

    # Linear fit with bounds
    ind2 = np.isfinite(rmse)
    def func(x,a,b):
        return a*x+b

    popt, pcov=curve_fit(func, dm[ind2]**2, rmse[ind2]**2, bounds=(0,np.inf))
    rmseFit = np.sqrt(func(dm**2,*popt))
    return rmseFit


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


def _noisemodel_binned(fraction, d, dm, w, sKnots):
    """Estimate residual spread without repeating bin statistics per pixel."""

    model_squared, residual_squared = _binned_noise_statistics(fraction, d, dm, w, sKnots)

    if model_squared.size < 2:
        raise ValueError("At least two populated geometry bins are required")

    def func(x, a, b):
        return a*x+b

    coefficients, _ = curve_fit(func, model_squared, residual_squared, bounds=(0, np.inf))

    return np.sqrt(func(dm**2, *coefficients))
