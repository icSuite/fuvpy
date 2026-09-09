#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:28:25 2022

@author: aohma
"""
import warnings

import numpy as np
import xarray as xr

from scipy.interpolate import BSpline
from scipy.optimize import curve_fit, lsq_linear
from scipy.stats import binned_statistic
from scipy.linalg import lstsq

from polplot import sdarngrid,bin_number
from fuvpy.utils import sh
from fuvpy.utils.sunlight import subsol


def _asymmetric_legacy_coordinate(legacy_coordinate, nightside_width):
    """Compress only the negative side of the legacy dayglow coordinate."""

    legacy_coordinate = np.asarray(legacy_coordinate)
    return np.where(legacy_coordinate >= 0, legacy_coordinate, nightside_width*np.tanh(legacy_coordinate/nightside_width))


def _asymmetric_legacy_knots(nightside_width):
    """Return compressed-nightside knots with strict-support endpoint padding."""

    knots = np.array([-nightside_width, -0.75*nightside_width, -0.25*nightside_width, 0, 0.25, 1.5, 3.5])
    # The transformed horizon can reach -width and ordinary pixels can equal
    # the upper anchor, while the existing fit mask uses strict inequalities.
    padding = 1e-6*(knots[-1] - knots[0])
    knots[0] -= padding
    knots[-1] += padding
    return knots


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


def _solve_monotonic_bs(weighted_design, weighted_data, damping, n_scp, n_tcp, solver):
    """Solve for spatially non-decreasing B-spline coefficients.

    BVLS is allowed at most 1000 iterations. SciPy's ``result.success`` is the
    convergence criterion; reaching the iteration limit is not interpreted
    separately.
    """

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

    lower = np.zeros(n_coefficients)
    lower[:n_tcp] = -np.inf

    try:
        result = lsq_linear(design, data, bounds=(lower, np.inf), method='bvls', max_iter=1000)
    except Exception as error:
        raise RuntimeError(f"Monotonic B-spline solver failed: {error}") from error
    if not result.success:
        raise RuntimeError(f"Monotonic B-spline fit failed (optimality={result.optimality:.3g}): {result.message}")

    return transform @ result.x

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
        Damping (Tikhonov regularization).
        The default is 0 (no damping).
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
    asymmetric_legacy : bool, optional
        If True, retain nonnegative `cos(SZA)/cos(DZA)` values and compress
        negative values as `nightside_width*tanh(x/nightside_width)`. Default
        is False.
    nightside_width : float, optional
        Positive finite scale for the asymmetric legacy nightside. Default is
        0.2.
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
        the legacy coordinate `cos(SZA)/cos(DZA)`. The temporal spline remains
        unchanged. Default is False.
    monotonic_solver : {'normal', 'full'}, optional
        `normal` first reduces the weighted image rows to a small normal
        system and is the operational default. `full` applies the same bounds
        to all weighted rows as a reference calculation.
    monotonic_fallback : bool, optional
        If True, a failed monotonic solve emits a warning and uses the ordinary
        solver for that IRLS iteration and all remaining iterations. Default
        is False, so monotonic solver failures remain explicit.
    diagnostics : bool, optional
        If True, add scalar `dgfit_iterations` and `dgfit_converged` fields.
        Default is False.
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
    noise_model = kwargs.pop('noise_model') if 'noise_model' in kwargs.keys() else 'legacy'
    spatial_weighting = kwargs.pop('spatial_weighting') if 'spatial_weighting' in kwargs.keys() else 'gaussian'
    asymmetric_legacy = bool(kwargs.pop('asymmetric_legacy')) if 'asymmetric_legacy' in kwargs.keys() else False
    nightside_width = kwargs.pop('nightside_width') if 'nightside_width' in kwargs.keys() else 0.2
    dgsigma_freeze_iteration = kwargs.pop('dgsigma_freeze_iteration') if 'dgsigma_freeze_iteration' in kwargs.keys() else 1
    input_variance = kwargs.pop('input_variance') if 'input_variance' in kwargs.keys() else None
    variance_weighted = bool(kwargs.pop('variance_weighted')) if 'variance_weighted' in kwargs.keys() else False
    monotonic = bool(kwargs.pop('monotonic')) if 'monotonic' in kwargs.keys() else False
    monotonic_solver = kwargs.pop('monotonic_solver') if 'monotonic_solver' in kwargs.keys() else 'normal'
    monotonic_fallback = bool(kwargs.pop('monotonic_fallback')) if 'monotonic_fallback' in kwargs.keys() else False
    diagnostics = bool(kwargs.pop('diagnostics')) if 'diagnostics' in kwargs.keys() else False
    inplace = bool(kwargs.pop('inplace')) if 'inplace' in kwargs.keys() else False

    if noise_model not in ['legacy', 'binned']:
        raise ValueError("noise_model must be 'legacy' or 'binned'")

    if spatial_weighting not in ['gaussian', 'bin_population']:
        raise ValueError("spatial_weighting must be 'gaussian' or 'bin_population'")

    try:
        nightside_width = float(nightside_width)
    except (TypeError, ValueError):
        raise ValueError("nightside_width must be a finite number greater than zero") from None
    if not np.isfinite(nightside_width) or nightside_width <= 0:
        raise ValueError("nightside_width must be a finite number greater than zero")

    if (isinstance(dgsigma_freeze_iteration, (bool, np.bool_)) or not isinstance(dgsigma_freeze_iteration, (int, np.integer)) or dgsigma_freeze_iteration < 0):
        raise ValueError("dgsigma_freeze_iteration must be a non-negative integer")

    if monotonic_solver not in ['normal', 'full']:
        raise ValueError("monotonic_solver must be 'normal' or 'full'")

    if not inplace: imgs = imgs.copy(deep=True)

    # Add temporal dimension if missing
    if len(imgs.sizes)==2: imgs = imgs.expand_dims('date')

    # Spatial knots and viewing angle correction
    if imgs['id'] in ['WIC','SI12','SI13','UVI']:
        fraction = np.cos(np.deg2rad(imgs['sza'].stack(z=('row','col')).values))/np.cos(np.deg2rad(imgs['dza'].stack(z=('row','col')).values))
        if asymmetric_legacy:
            fraction = _asymmetric_legacy_coordinate(fraction, nightside_width)
        if sKnots is None:
            sKnots = _asymmetric_legacy_knots(nightside_width) if asymmetric_legacy else [-3.5,-0.25,0,0.25,1.5,3.5]
    elif imgs['id'] == 'VIS':
        if asymmetric_legacy:
            raise ValueError("asymmetric_legacy requires the legacy FUV coordinate and is not available for VIS")
        fraction = np.cos(np.deg2rad(imgs['sza'].stack(z=('row','col')).values))
        if sKnots is None: sKnots= [-1,-0.1,0,0.1,0.4,1]
        if spatial_weighting == 'bin_population':
            raise ValueError("bin_population spatial weighting requires the legacy FUV coordinate")
        if monotonic:
            raise ValueError("monotonic fitting requires the legacy FUV coordinate")
    sKnots = np.r_[np.repeat(sKnots[0],sOrder),sKnots, np.repeat(sKnots[-1],sOrder)]

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
    n_scp = len(sKnots)-sOrder-1

    # Temporal design matix
    M = BSpline(tKnots, np.eye(n_tcp), tOrder)(time)

    # Spatial design matrix
    G_g=[]
    G_s=[]
    for i in range(n_t):
        G= BSpline(sKnots, np.eye(n_scp), sOrder)(fraction[i,:]) # Spatial design matirx
        G_t = np.zeros((n_s, n_scp*n_tcp))
        for j in range(n_tcp):
            G_t[:, np.arange(j, n_scp*n_tcp, n_tcp)] = G*M[i, j]

        G_g.append(G)
        G_s.append(G_t)
    G_s=np.vstack(G_s)

    # Index of data to use in model
    ind = (imgs['sza'].stack(z=('row','col')).values >= 0) & (np.isfinite(imgs[inImg].stack(z=('row','col')).values)) & (fraction>sKnots[0]) & (fraction<sKnots[-1])

    # Spatial sampling weights
    if spatial_weighting == 'gaussian':
        ws = 1/np.exp(-0.5*(fraction/1.5)**2)
    else:
        spatial_bins = np.arange(sKnots[0], sKnots[-1] + 0.1, 0.1)
        ws = _bin_population_weights(fraction, ind, spatial_bins)

    # Make everything flat
    d_s = imgs[inImg].stack(z=('row','col')).values.flatten()
    ind = ind.flatten()
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

    # Damping (zeroth-order Tikhonov regularization)
    damping = dampingVal*np.ones(G_s.shape[1])
    R = np.diag(damping)

    # Iterativaly solve the inverse problem
    diff = 1e10
    iteration = 0
    m = None
    use_monotonic = monotonic
    solver_flag = 1 if monotonic else 0
    while (diff>stop)&(iteration < 100):
        # Keep robust, spatial-sampling, and measurement-precision roles
        # visible before forming the weighted least-squares system.
        fit_w = w[ind] * ws[ind] * measurement_w[ind]
        weighted_design = G_s[ind, :] * fit_w[:, None]
        weighted_data = d_s[ind] * fit_w
        if use_monotonic:
            try:
                m_s = _solve_monotonic_bs(weighted_design, weighted_data, R, n_scp, n_tcp, monotonic_solver)
            except RuntimeError as error:
                if not monotonic_fallback:
                    raise
                warnings.warn(f"Monotonic B-spline solve failed; using the ordinary solver for this and all remaining IRLS iterations: {error}", RuntimeWarning, stacklevel=2)
                use_monotonic = False
                solver_flag = 2
        if not use_monotonic:
            m_s = lstsq(weighted_design.T@weighted_design + R, weighted_design.T@weighted_data, lapack_driver='gelsy')[0]
        mNew    = M@m_s.reshape((n_scp, n_tcp)).T
        dm=[]
        for i, tt in enumerate(time):
            dm.append(G_g[i]@mNew[i, :])

        dm=np.array(dm).squeeze()
        residuals = (d_s[ind] - dm.flatten()[ind])

        if iteration <= dgsigma_freeze_iteration:
            sigma = np.full_like(d_s,np.nan)
            if noise_model == 'legacy':
                sigma[ind] = _noisemodel(fraction.flatten()[ind], d_s[ind], dm.flatten()[ind], w[ind]*ws[ind], sKnots)
            else:
                sigma[ind] = _noisemodel_binned(fraction.flatten()[ind], d_s[ind], dm.flatten()[ind], w[ind]*ws[ind], sKnots)
        w[ind] = (1 - np.minimum(1,(residuals/(tukeyVal*sigma[ind]))**2))**2

        if m is not None: diff = np.sqrt(np.mean((mNew-m)**2))/(1+np.sqrt(np.mean(mNew**2)))
        m = mNew
        iteration += 1

    if iteration == 100: print('Warning: Model did not converge. Truncated after 100 iterations')

    normM = np.sqrt(np.average(m_s**2))
    normR = np.sqrt(np.average(residuals**2,weights=w[ind]))

    # Add dayglow model and corrected image to the Dataset
    imgs['dgmodel'] = (['date','row','col'],(dm).reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['dgimg'] = imgs[inImg]-imgs['dgmodel']
    imgs['dgweight'] = (['date','row','col'],(w).reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['dgsigma'] = (['date','row','col'],(sigma).reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['dgnorm_residual'] = normR
    imgs['dgnorm_model'] = normM
    imgs['dgsolver_flag'] = solver_flag
    if diagnostics:
        imgs['dgfit_iterations'] = iteration
        imgs['dgfit_converged'] = iteration < 100 and np.isfinite(diff) and diff <= stop

    # Remove pixels outside model scope
    ind = (imgs.sza>=0)
    imgs['dgmodel'] = xr.where(~ind,np.nan,imgs['dgmodel'])
    imgs['dgimg'] = xr.where(~ind,np.nan,imgs['dgimg'])
    imgs['dgweight'] = xr.where(~ind,np.nan,imgs['dgweight'])
    imgs['dgsigma'] = xr.where(~ind,np.nan,imgs['dgsigma'])

    # Add attributes
    imgs['dgmodel'].attrs = {'long_name': 'BS model'}
    imgs['dgimg'].attrs = {'long_name': 'BS corrected image'}
    imgs['dgweight'].attrs = {'long_name': 'BS model weights'}
    imgs['dgsigma'].attrs = {'long_name': 'BS model spread'}
    imgs['dgnorm_residual'].attrs = {'long_name': 'BS model residual norm'}
    imgs['dgnorm_model'].attrs = {'long_name': 'BS model solution norm'}
    imgs['dgsolver_flag'].attrs = {
        'long_name': 'BS spatial solver outcome',
        'flag_values': np.array([0, 1, 2], dtype=np.int8),
        'flag_meanings': 'standard_requested monotonic_succeeded monotonic_failed_standard_fallback'
    }
    if diagnostics:
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


def backgroundmodel_SH(imgs,Nsh,Msh,**kwargs):
    '''
    Function to model the FUV residual background and subtract it from the input image

    Parameters
    ----------
    imgs : xarray.Dataset
        Dataset with the FUV images
    Nsh : int
        Order of the SH
    Msh : int
        Degree of the SH
    dampingVal : TYPE, optional
        Damping to reduce the influence of the time-dependent part of the model.
        The default is 0 (no damping).
    tukeyVal : float, optional
        Determines to what degree outliers are down-weighted.
        Iterative reweights is (1-(residuals/(tukeyVal*rmse))^2)^2
        Larger tukeyVal means less down-weight
        Default is 5
    stop : float, optional
        When to stop the iteration. The default is 0.001.
    n_tKnots : int, optional
        Number of temporal knots, equally spaced between the endpoints. The default is 2 (only knots at endpoints)
    tOrder: int, optional
        Order of the temporal spline fit. The default is 2.
    inplace : bool, optional
        If True, update the Dataset in place. Default is False (return new Dataset)
    Returns
    -------
    imgs : xarray.Dataset, if inplace is False
        A copy of the input Dataset with new data_vars:
            - imgs['shmodel'] is the sh model
            - imgs['shimg'] is the sh-corrected image
            - imgs['shweight'] is the weights after the final iteration
            - imgs['shnorm_residual'] is the residual norm
            - imgs['shnorm_model'] is the model norm
    '''

    # Set keyword arguments to input or default values    
    dampingVal = kwargs.pop('dampingVal') if 'dampingVal' in kwargs.keys() else 0
    tukeyVal = kwargs.pop('tukeyVal') if 'tukeyVal' in kwargs.keys() else 5
    stop = kwargs.pop('stop') if 'stop' in kwargs.keys() else 1e-3
    n_tKnots = kwargs.pop('n_tKnots') if 'n_tKnots' in kwargs.keys() else 2
    tOrder = kwargs.pop('tOrder') if 'tOrder' in kwargs.keys() else 2
    inplace = bool(kwargs.pop('inplace')) if 'inplace' in kwargs.keys() else False

    if not inplace: imgs = imgs.copy(deep=True)

    # Add temporal dimension if missing
    if len(imgs.sizes)==2: imgs = imgs.expand_dims('date')

    time=(imgs['date'].values-imgs['date'].values[0])/ np.timedelta64(1, 'm')

    glat = imgs['glat'].stack(z=('row','col')).values
    d = imgs['dgimg'].stack(z=('row','col')).values

    # Normalize on dayglow noise
    d = d/imgs['dgsigma'].stack(z=('row','col')).values

    sslat, sslon = map(np.ravel, subsol(imgs['date'].values))
    phi = np.deg2rad((imgs['glon'].stack(z=('row','col')).values - sslon[:,None] + 180) % 360 - 180)

    n_t = glat.shape[0]
    n_s = glat.shape[1]

    # Temporal knots
    tKnots = np.linspace(time[0], time[-1], n_tKnots)
    tKnots = np.r_[np.repeat(tKnots[0],tOrder),tKnots, np.repeat(tKnots[-1],tOrder)]

    # Number of control points
    n_cp = len(tKnots)-tOrder-1

    # Temporal design matix
    M = BSpline(tKnots, np.eye(n_cp), tOrder)(time)

    # Iterative (few iterations)
    skeys = sh.SHkeys(Nsh, Msh).Mge(1).MleN().NminusMeven()
    ckeys = sh.SHkeys(Nsh, Msh).MleN().NminusMeven()

    G_g=[]
    G_s=[]
    for i in range(n_t):
        # calculate Legendre functions at glat:
        P, dP = sh.get_legendre(Nsh, Msh, 90 - glat[i,:])
        Pc = np.hstack([P[key] for key in ckeys])
        Ps = np.hstack([P[key] for key in skeys])

        Gcos = Pc * np.cos(phi[i,:].reshape((-1, 1)) * ckeys.m)
        Gsin = Ps * np.sin(phi[i,:].reshape((-1, 1)) * skeys.m)

        G = np.hstack((Gcos, Gsin))
        G = G/imgs['dgsigma'].stack(z=('row','col')).values[i,:][:,None]
        n_G = np.shape(G)[1]

        G_t = np.zeros((n_s, n_G*n_cp))

        for j in range(n_cp):
            G_t[:, np.arange(j, n_G*n_cp, n_cp)] = G*M[i, j]

        G_g.append(G)
        G_s.append(G_t)

    G_s = np.array(G_s)
    G_s = G_s.reshape(-1,G_s.shape[2])

    # Data
    ind = (np.isfinite(d))&(glat>0)

    # Spatial weights
    grid,mltres=sdarngrid(5,5,np.nanmin(glat)//5*5) # Equal area grid
    ws = np.full(glat.shape,np.nan)
    for i in range(len(glat)):
        gbin = bin_number(grid,glat[i,ind[i]],imgs['glon'].stack(z=('row','col')).values[i,ind[i]]/15)
        count=np.full(grid.shape[1],0)
        un,count_un=np.unique(gbin,return_counts=True)
        count[un]=count_un
        ws[i,ind[i]]=1/count[gbin]

    # Make everything flat
    d_s = d.flatten()
    ind = ind.flatten()
    ws = ws.flatten()
    w = imgs['dgweight'].stack(z=('row','col')).values.flatten()

    # Damping
    damping = dampingVal*np.ones(G_s.shape[1])
    R = np.diag(damping)

    diff = 1e10*stop
    iteration = 0
    m = None
    while (diff>stop)&(iteration < 100):
        # Solve for spline amplitudes
        m_s = np.linalg.lstsq((G_s[ind,:]*w[ind,None]*ws[ind,None]).T@(G_s[ind,:]*w[ind,None]*ws[ind,None])+R,(G_s[ind,:]*w[ind,None]*ws[ind,None]).T@(d_s[ind]*w[ind]*ws[ind]),rcond=None)[0]

        # Retrieve B-spline smooth model paramters (coarse)
        mNew    = M@m_s.reshape((n_G, n_cp)).T
        dm=[]
        for i, tt in enumerate(time):
            dm.append(G_g[i]@mNew[i, :])

        dm=np.array(dm).squeeze()
        residuals = dm.flatten()[ind] - d.flatten()[ind]

        w[ind] = (1 - np.minimum(1,(residuals/(tukeyVal))**2))**2

        if m is not None:
            diff = np.sqrt(np.mean( (mNew-m)**2))/(1+np.sqrt(np.mean(mNew**2)))
        m = mNew
        iteration += 1

    if iteration == 100: print('Warning: Model did not converge. Truncated after 100 iterations')
    normM = np.sqrt(np.average(m_s**2))
    normR = np.sqrt(np.average(residuals**2,weights=w[ind]))

    imgs['shmodel'] = (['date','row','col'],(dm*imgs['dgsigma'].stack(z=('row','col')).values).reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['shimg'] = imgs['dgimg']-imgs['shmodel']
    imgs['shweight'] = (['date','row','col'],(w).reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['shnorm_residual'] = normR
    imgs['shnorm_model'] = normM

    # Remove pixels outside model scope
    ind = (imgs.glat >= 0) & np.isfinite(imgs['dgsigma'])
    imgs['shmodel'] = xr.where(~ind,np.nan,imgs['shmodel'])
    imgs['shimg'] = xr.where(~ind,np.nan,imgs['shimg'])
    imgs['shweight'] = xr.where(~ind,np.nan,imgs['shweight'])

    # Add attributes
    imgs['shmodel'].attrs = {'long_name': 'SH model'}
    imgs['shimg'].attrs = {'long_name': 'SH corrected image'}
    imgs['shweight'].attrs = {'long_name': 'SH model weights'}
    imgs['shnorm_residual'].attrs = {'long_name': 'SH model residual norm'}
    imgs['shnorm_model'].attrs = {'long_name': 'SH model solution norm'}

    if not inplace:
        return imgs


def bin_number2(grid, mlat, mlt):
    """
    Ultra-fast routine to determine the bin number of mlat and mlt
    NB HARD CODED VERSION TO USE DATA FROM BOTH HEMISPHERES

    Parameters
    ----------
    grid : 2 x N array
        Array containing the edges of the mlat-mlt bins.
        First row is mlat and second row mlt.
        grid can be constructed using either equal_area_grid or sdarngrid
    mlat : array
        Array with the mlats to bin
    mlt : array
        Array with the mlts to bin

    Returns
    -------
    bin_n : array
        Array with the bin number for each mlat-mlt pair.
        Locations outside the defined grid are set to -1 in the returned array.


    SMH 2021/04/15
    Modified by JPR 2021/10/20
    2022-05-02: JPR added nan handling
    2022-05-05: AO added better handling of values outside grid (including nans).
    """

    llat = np.unique(grid[0]) # latitude circles
    assert np.allclose(np.sort(llat) - llat, 0) # should be in sorted order automatically. If not, the algorithm will not work
    dlat = np.diff(llat)[0] # latitude step
    latbins = np.hstack(( llat, llat[-1] + dlat )) # make latitude bin edges

    bin_n = -np.ones_like(mlat).astype(int) # initiate binnumber
    ii = (np.isfinite(mlat))&(np.isfinite(mlt))&(mlat>=latbins[0])&(mlat<=latbins[-1]) # index of real values inside the grid
    mlat,mlt = mlat[ii],mlt[ii] # Reduce to only real values inside grid

    latbin_n = np.digitize(mlat, latbins) - 1 # find the latitude index for each data point

    # number of longitude bins in each latitude ring:
    nlons = np.array([len(np.unique(grid[1][grid[0] == lat])) for lat in llat])

    # normalize all longitude bins to the equatorward ring:
    _mlt = mlt * nlons[latbin_n] / nlons[17]

    # make longitude bin edges for the equatorward ring:
    llon = np.unique(grid[1][grid[0] == llat[17]])
    dlon = np.diff(llon)[0]
    lonbins = np.hstack((llon, llon[-1] + dlon)) # make longitude bin edges
    lonbin_n = np.digitize(_mlt, lonbins) - 1 # find the longitude bin

    # map from 2D bin numbers to 1D by adding the number of bins in each row equatorward:
    bin_n[ii] = lonbin_n + np.cumsum(np.hstack((0, nlons)))[latbin_n]

    return bin_n

def backgroundmodel_SBS(imgs,**kwargs):
    '''
    Model background emissions using Spherical B-spline (SBS) modeling

    Parameters
    ----------
    imgs : xarray.Dataset
        Dataset with the FUV images, imported by readFUVimage()
    inImg : str, optional
        Name of the input image to be used in the model. The default is 'img'.
    latOrder : int, optional
        Order of the latitudinal B-splines. The default is 3.
    lonOrder : int, optional
        Order of the longitudinal B-splines. The default is 3.
    tOrder : int, optional
        Order of the temporal B-splines. The default is 2.
    latKnots : array like, optional
        Location of the latitudinal knots. The default is None (default is used).
    lonKnots : array like, optional
        Location of the longitudinal knots. The default is None (default is used).
    n_tKnots : int, optional
        Number of temporal knots, equally spaced between the endpoints. The default is 2 (only knots at endpoints)

    dampingVal : float, optional
        Damping (Tikhonov regularization).
        The default is 0 (no damping).
    tukeyVal : float, optional
        Determines to what degree outliers are down-weighted.
        Iterative reweights is (1-(residuals/(tukeyVal*rmse))^2)^2
        Larger tukeyVal means less down-weight
        Default is 5
    stop : float, optional
        When to stop the iteration. The default is 0.001.
    height : float, optional
        Assumed height of the emissions in km. Default is 130
    dzacorr : float, optional
        Absorption constant for viewning angle correction. Default is 0.15
    inplace : bool, optional
        If True, update the Dataset in place. Default is False (return a new Dataset).


    Returns
    -------
    imgs : xarray.Dataset, if inplace is False
        A copy of the input Dataset with new data_vars:
            - imgs['dgmodel'] is the B-spline based dayglow model
            - imgs['dgimg'] is the dayglow-corrected image (dayglow subtracked from the input image)
            - imgs['dgweight'] is the weights after the final iteration
            - imgs['dgnorm_residual'] is the residual norm
            - imgs['dgnorm_model'] is the model norm
            - imgs['background'] is the identified background level
    '''

    # Set keyword arguments to input or default values    
    inImg = kwargs.pop('inImg') if 'inImg' in kwargs.keys() else 'img'
    latOrder = kwargs.pop('latOrder') if 'latOrder' in kwargs.keys() else 3
    lonOrder = kwargs.pop('lonOrder') if 'lonOrder' in kwargs.keys() else 3
    tOrder = kwargs.pop('tOrder') if 'tOrder' in kwargs.keys() else 2
    latKnots = kwargs.pop('latKnots') if 'latKnots' in kwargs.keys() else [-90,-10,0,10,50,90]
    lonKnots = kwargs.pop('lonKnots') if 'lonKnots' in kwargs.keys() else np.array([0,90,180,270])
    n_tKnots = kwargs.pop('n_tKnots') if 'n_tKnots' in kwargs.keys() else 2

    dampingVal = kwargs.pop('dampingVal') if 'dampingVal' in kwargs.keys() else 0
    tukeyVal = kwargs.pop('tukeyVal') if 'tukeyVal' in kwargs.keys() else 5
    stop = kwargs.pop('stop') if 'stop' in kwargs.keys() else 1e-3
    height = kwargs.pop('height') if 'height' in kwargs.keys() else 130
    dzacorr = kwargs.pop('dzacorr') if 'dzacorr' in kwargs.keys() else 0.15
    
    inplace = bool(kwargs.pop('inplace')) if 'inplace' in kwargs.keys() else False

    if not inplace: imgs = imgs.copy(deep=True)

    sslat, sslon = map(np.ravel, subsol(imgs['date'].values))


    r = 6371 + height
    theta = np.deg2rad(90-imgs.glat.values)
    phi = np.deg2rad((imgs['glon'].values - sslon[:,None,None] + 180) % 360 - 180)

    xGEO = r * np.sin(theta) * np.cos(phi)
    yGEO = r * np.sin(theta) * np.sin(phi)
    zGEO = r * np.cos(theta)

    slat = np.full_like(theta,np.nan)
    slon = np.full_like(theta,np.nan)
    for i in range(len(imgs.date)):
        conv = np.array([[-np.sin(np.deg2rad(sslat[i])),0,np.cos(np.deg2rad(sslat[i]))],
                         [0,-1,0],
                         [np.cos(np.deg2rad(sslat[i])),0,np.sin(np.deg2rad(sslat[i]))]])
        xyz = conv @ np.vstack((xGEO[i,:,:].flatten(),yGEO[i,:,:].flatten(),zGEO[i,:,:].flatten()))

        slat[i,:,:] = 90-np.rad2deg(np.arccos(xyz[2,:]/r)).reshape(xGEO[i,:,:].shape)
        np.seterr(invalid='ignore', divide='ignore')
        slon[i,:,:] = np.rad2deg(((np.arctan2(xyz[1,:], xyz[0,:])*180/np.pi) % 360)/180*np.pi).reshape(xGEO[i,:,:].shape)

    imgs['slat'] = (['date','row','col'],slat)
    imgs['slon'] = (['date','row','col'],slon)

    # Viewving angle correction
    background=np.nanmedian(imgs[inImg].values[(imgs['sza'].values>100)|(np.isnan(imgs['sza'].values))])
    imgs['cimg'] = (imgs[inImg]-background)*np.cos(np.deg2rad(imgs['dza']))/np.exp(dzacorr*(1. - 1/np.cos(np.deg2rad(imgs['dza']))))

    # Coordinates and data
    time=(imgs['date'].to_series().values-imgs['date'].values[0])/ np.timedelta64(1, 'm')
    slat = imgs['slat'].stack(z=('row','col')).values
    slon = imgs['slon'].stack(z=('row','col')).values
    d = imgs['cimg'].stack(z=('row','col')).values

    ind = (imgs['sza'].stack(z=('row','col')).values >= 0) & (np.isfinite(imgs['cimg'].stack(z=('row','col')).values))

    # temporal and spatial size
    n_t = slat.shape[0]

    # colat knots
    latKnots = np.r_[np.repeat(latKnots[0],latOrder),latKnots, np.repeat(latKnots[-1],latOrder)]
    n_latcp = len(latKnots)-latOrder-1 # Number of control points

    # lon knots full
    lonKnots = np.r_[lonKnots-360,lonKnots,lonKnots+360]
    lonKnots = np.r_[np.repeat(lonKnots[0],lonOrder),lonKnots, np.repeat(lonKnots[-1],lonOrder)]
    n_loncp = len(lonKnots)-lonOrder-1 # Number of control points

    # Temporal knots
    if n_t==1:
        tKnots = np.linspace(time[0], time[-1]+1, n_tKnots)
        tOrder = 0
    else:
        tKnots = np.linspace(time[0], time[-1], n_tKnots)
        tKnots = np.r_[np.repeat(tKnots[0],tOrder),tKnots, np.repeat(tKnots[-1],tOrder)]
    tKnots = np.r_[np.repeat(tKnots[0],tOrder),tKnots, np.repeat(tKnots[-1],tOrder)]
    n_tcp = len(tKnots)-tOrder-1 # Number of control points

    # Temporal design matix
    M = BSpline(tKnots, np.eye(n_tcp), tOrder)(time)

    G_g=[]
    G_s=[]
    for t in range(n_t):
        # Colat design matix
        n_s = np.sum(ind[t,:]) # Spatial data points at time t
        Glat = BSpline(latKnots, np.eye(n_latcp), latOrder)(slat[t,ind[t,:]])
        Gtemp = BSpline(lonKnots, np.eye(n_loncp), lonOrder)(slon[t,ind[t,:]])
        Glon = Gtemp[:,4:8].copy()+Gtemp[:,8:12].copy()

        G = np.tile(Glon,(1,n_latcp))*np.repeat(Glat,4,axis=1)
        G_t = np.tile(G,(1,n_tcp))*np.repeat(M[t,:],n_latcp*4)[None,:]

        G = np.tile(Glon,(1,n_latcp))*np.repeat(Glat,4,axis=1)
        G_t = np.repeat(G,n_tcp,axis=1)*np.tile(M[t,:],n_latcp*4)[None,:]

        G_g.append(G)
        G_s.append(G_t)
    G_s=np.vstack(G_s)

    # Spatial weights
    grid,mltres=sdarngrid(5,5,-90) # Equal area grid
    ws = np.full(slat.shape,np.nan)
    for i in range(len(slat)):
        sbin = bin_number2(grid,slat[i,ind[i]],slon[i,ind[i]]/15)
        count=np.full(grid.shape[1],0)
        un,count_un=np.unique(sbin,return_counts=True)
        count[un]=count_un
        ws[i,ind[i]]=1/count[sbin]

    # Make everything flat
    d_s = d[ind]
    ws = ws[ind]
    ws[:]=1 # TEMP NO SPATIAL WEIGHTING
    w = np.ones_like(d_s)

    # Damping (zeroth-order Tikhonov regularization)
    damping = dampingVal**np.ones(G_s.shape[1])
    R = np.diag(damping)

    # 1- order Tikhonov in the longitudinal direction 
    L = np.array([[-1,1,0,0],[0,-1,1,0],[0,0,-1,1],[1,0,0,-1]])
    LTL = np.zeros((4*n_tcp*n_latcp,4*n_tcp*n_latcp))
    for i in range(n_latcp):
        for t in range(n_tcp):
            LTL[t+(i*4*n_tcp):t+4*n_tcp+(i*4*n_tcp):n_tcp,t+(i*4*n_tcp):t+4*n_tcp+(i*4*n_tcp):n_tcp] = L.T@L
    R = dampingVal*LTL

    # Iterativaly solve the inverse problem
    diff = 1e10
    iteration = 0
    m = None

    while (diff>stop)&(iteration < 100):
        m_s = np.linalg.lstsq((G_s*w[:,None]*ws[:,None]).T@(G_s*w[:,None]*ws[:,None])+R,(G_s*w[:,None]*ws[:,None]).T@(d_s*w*ws),rcond=None)[0]

        mNew    = M@m_s.reshape((n_latcp*4, n_tcp)).T
        dm=[]
        for i, tt in enumerate(time):
            dm.append((G_g[i]@mNew[i, :]).flatten())
        dm = np.hstack(dm)
        residuals = (d_s - dm)

        rmse = np.sqrt(np.average(residuals**2)     )
        w = (1 - np.minimum(1,(residuals/(tukeyVal*rmse))**2))**2

        if m is not None: diff = np.sqrt(np.mean((mNew-m)**2))/(1+np.sqrt(np.mean(mNew**2)))

        m = mNew
        iteration += 1

    normM = np.sqrt(np.average(m_s**2))
    normR = np.sqrt(np.average(residuals**2,weights=w))

    # Add dayglow model and corrected image to the Dataset
    model = np.full_like(d,np.nan)
    model[ind]= dm
    model = model.reshape((n_t,len(imgs.row),len(imgs.col))) * np.exp(dzacorr*(1. - 1/np.cos(np.deg2rad(imgs['dza'])))) / np.cos(np.deg2rad(imgs['dza'])) + background

    weights = np.full_like(d,np.nan)
    weights[ind]=w
    imgs['dgmodel'] = model
    imgs['dgimg'] = imgs['img']-imgs['dgmodel']
    imgs['dgweight'] = (['date','row','col'],(weights).reshape((n_t,len(imgs.row),len(imgs.col))))
    imgs['dgnorm_residual'] = normR
    imgs['dgnorm_model'] = normM
    imgs['background'] = background

    # Remove pixels outside model scope
    ind = (imgs.sza>=0)
    imgs['dgmodel'] = xr.where(~ind,np.nan,imgs['dgmodel'])
    imgs['dgimg'] = xr.where(~ind,np.nan,imgs['dgimg'])
    imgs['dgweight'] = xr.where(~ind,np.nan,imgs['dgweight'])

    # Add attributes
    imgs['dgmodel'].attrs = {'long_name': 'BS model'}
    imgs['dgimg'].attrs = {'long_name': 'BS corrected image'}
    imgs['dgweight'].attrs = {'long_name': 'BS model weights'}
    imgs['dgnorm_residual'].attrs = {'long_name': 'BS model residual norm'}
    imgs['dgnorm_model'].attrs = {'long_name': 'BS model solution norm'}
    imgs['background'].attrs = {'long_name': 'Background level'}
    
    # Return the new DataSet if not inplace 
    if not inplace:
        return imgs
