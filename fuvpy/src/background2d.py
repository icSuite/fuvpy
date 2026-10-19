"""Monotonic x/solar-view-sine backgrounds, temporal then independent."""

from collections.abc import Mapping

import numpy as np
from scipy.interpolate import BSpline
import xarray as xr

from fuvpy.src import background as bg


__all__ = ['backgroundmodel_BS2D']


def _spline(breaks, name):
    breaks = np.asarray(breaks, dtype=float)
    if breaks.ndim != 1 or len(breaks) < 2 or np.any(~np.isfinite(breaks)) or np.any(np.diff(breaks) <= 0):
        raise ValueError(f'{name} must be finite, strictly increasing break points')
    knots = np.r_[np.repeat(breaks[0], 3), breaks, np.repeat(breaks[-1], 3)]
    return BSpline(knots, np.eye(len(knots)-4), 3, extrapolate=False)


def _reference_ridge(images, temporal, directional, damping, variance_weighted):
    """Keep the validated numerical ridge scale of the legacy reference basis.

    This design is used only to set a scalar penalty, never as a fitted term.
    Its base/directional time targets are120/60 minutes at default settings.
    """
    counts = images.img.values.reshape(images.sizes['date'], -1)
    x = (np.cos(np.deg2rad(images.sza.values))/np.cos(np.deg2rad(images.dza.values))).reshape(counts.shape)
    q = images.solar_view_sine_factor.values.reshape(counts.shape)
    taper = bg._spencer_taper_weights(images.sza.values.reshape(counts.shape))
    valid = np.isfinite(counts) & (images.sza.values.reshape(counts.shape) >= 0) & (x > -3.5) & (x < 3.5) & np.isfinite(q) & np.isfinite(taper)
    population = bg._bin_population_weights(x, valid, np.arange(-3.5, 3.5+.1, .1))
    precision = np.ones_like(counts, dtype=counts.dtype if counts.dtype.kind == 'f' else float)
    if variance_weighted:
        variance = images.img_variance.values.reshape(counts.shape)
        if np.any(~np.isfinite(variance[valid])) or np.any(variance[valid] <= 0):
            raise ValueError('img_variance must be positive and finite at fitted pixels')
        precision[valid] = 1/np.sqrt(variance[valid])
    splines = [_spline([-3.5, -.25, 0, .25, 1.5, 3.5], 'reference x'), _spline([-3.5, 0, 1.5, 3.5], 'reference directional x')]
    diagonal = [np.zeros((s.c.shape[0], t.shape[1])) for s, t in zip(splines, [temporal, directional])]
    for frame in range(len(counts)):
        use = valid[frame]
        base, signed = [s(x[frame, use]) for s in splines]
        signed *= (q[frame, use]*taper[frame, use])[:, None]
        row_weight = (population[frame, use]*precision[frame, use])[:, None]
        for part, design, time_basis in zip(diagonal, [base, signed], [temporal, directional]):
            part += np.sum((design*row_weight)**2, axis=0)[:, None]*time_basis[frame]**2
    normal_diagonal = np.concatenate([part.ravel() for part in diagonal])
    return float(damping*bg._normal_matrix_scale(np.diag(normal_diagonal)))


def _fit_sequence(images, x_spline, sine_spline, temporal, ridge, freeze, stop, tukey, maximum, variance_weighted, sigma_coefficients=None):
    """Stream spatial normal systems; constrain every sine/time x sequence."""
    dtype = images.img.dtype if images.img.dtype.kind == 'f' else np.dtype(float)
    counts = images.img.values.astype(dtype, copy=False).reshape(images.sizes['date'], -1)
    sza = images.sza.values.reshape(counts.shape)
    x = np.cos(np.deg2rad(sza))/np.cos(np.deg2rad(images.dza.values.reshape(counts.shape)))
    q = images.solar_view_sine_factor.values.reshape(counts.shape)
    lower, upper = x_spline.t[0], x_spline.t[-1]
    valid = np.isfinite(counts) & (sza >= 0) & (x > lower) & (x < upper) & np.isfinite(q) & (np.abs(q) <= 1)
    if np.any(~valid.any(axis=1)):
        raise ValueError('Every selected frame must contain supported x/sine pixels')
    indices = [np.flatnonzero(mask) for mask in valid]
    population = bg._bin_population_weights(x, valid, np.arange(lower, upper+.1, .1))
    precision = np.ones_like(counts)
    if variance_weighted:
        variance = images.img_variance.values.reshape(counts.shape)
        if np.any(~np.isfinite(variance[valid])) or np.any(variance[valid] <= 0):
            raise ValueError('img_variance must be positive and finite at fitted pixels')
        precision[valid] = 1/np.sqrt(variance[valid])
    n_x, n_sine, n_time = x_spline.c.shape[0], sine_spline.c.shape[0], temporal.shape[1]
    spatial_count = n_x*n_sine
    basis_scale = np.sqrt(n_sine)

    def design(frame):
        use = indices[frame]
        return basis_scale*(x_spline(x[frame, use])[:, :, None]*sine_spline(q[frame, use])[:, None, :]).reshape(len(use), spatial_count)

    # Fresh robust weights in both stages. Only the noise relation is shared.
    estimate_sigma = sigma_coefficients is None
    sigma = np.full_like(counts, np.nan)
    weight = np.ones_like(counts)
    damping = ridge*np.eye(spatial_count*n_time)
    previous, history = None, []
    for iteration in range(1, maximum+1):
        normal = np.zeros_like(damping)
        rhs = np.zeros(damping.shape[0])
        for frame, use in enumerate(indices):
            spatial = design(frame)
            row_weight = np.sqrt(weight[frame, use])*population[frame, use]*precision[frame, use]
            weighted = spatial*row_weight[:, None]
            local_normal = weighted.T@weighted
            local_rhs = weighted.T@(counts[frame, use]*row_weight)
            active = np.flatnonzero(temporal[frame])
            for a in active:
                rhs[a::n_time] += temporal[frame, a]*local_rhs
                for b in active:
                    normal[a::n_time, b::n_time] += temporal[frame, a]*temporal[frame, b]*local_normal
        coefficients = bg._solve_monotonic_bs(normal, rhs, damping, n_x, n_sine*n_time, 'normal_system')
        amplitudes = temporal@coefficients.reshape(spatial_count, n_time).T
        model = np.full(counts.shape, np.nan)
        for frame, use in enumerate(indices):
            model[frame, use] = design(frame)@amplitudes[frame]
        if estimate_sigma and iteration <= freeze:
            estimated, sigma_coefficients = bg._noisemodel_binned(x[valid], counts[valid], model[valid], (weight*population**2)[valid], x_spline.t, return_coefficients=True)
            sigma[valid] = estimated
            history.append(np.asarray(sigma_coefficients).copy())
        elif not estimate_sigma:
            sigma[valid] = np.sqrt(sigma_coefficients[0]*model[valid]**2 + sigma_coefficients[1])
        if np.any(~np.isfinite(sigma[valid])) or np.any(sigma[valid] <= 0):
            raise ValueError('The residual sigma model must give positive finite scales')
        residual = counts-model
        weight[valid] = (1-np.minimum(1, (residual[valid]/(tukey*sigma[valid]))**2))**2
        if not np.any(weight[valid] > 0):
            raise RuntimeError('2D background rejected all supported pixels')
        change = np.inf if previous is None else np.sqrt(np.mean((amplitudes-previous)**2))/(1+np.sqrt(np.mean(amplitudes**2)))
        previous = amplitudes.copy()
        if change <= stop and (not estimate_sigma or iteration > freeze):
            break
    else:
        raise RuntimeError(f'2D background did not converge after {maximum} fits (relative change {change:.6g})')
    shape = images.img.shape
    return {
        'model': model.reshape(shape), 'img': residual.reshape(shape), 'sigma': sigma.reshape(shape),
        'weight': np.where(valid, weight, np.nan).reshape(shape),
        'coefficients': coefficients.reshape(n_x, n_sine, n_time), 'sigma_coefficients': np.asarray(sigma_coefficients),
        'sigma_history': np.asarray(history), 'iterations': iteration, 'change': float(change)
    }


def backgroundmodel_BS2D(
    imgs, *, sKnots=None, sine_knots=None, temporal_knot_spacing=120, tOrder=2,
    dgsigma_freeze_iteration=2, stop=.01, dampingVal=.01, tukeyVal=5,
    variance_weighted=True, minimum_frame_quality=2, max_iterations=100,
    pointing_correction=False, pointing_options=None
):
    """Fit temporal and independent monotonic x/sine tensor backgrounds.

    ``x=cos(SZA)/cos(DZA)``; the second coordinate is exactly the stored
    ``solar_view_sine_factor``, including its illumination/view amplitudes.
    Both axes are cubic. Default x breaks give13 controls; [-1,1] gives four
    sine controls. There is no additional asymmetry or Spencer-tapered term.

    Estimate pooled residual sigma after each fit through
    ``dgsigma_freeze_iteration`` (default2), then freeze that pixel field and
    converge the temporal fit. Automatically fit each frame independently,
    sharing only the final ``sigma²=a*background²+b`` relation, evaluated at
    its own background each iteration. Tukey weights always start fresh.

    Temporal fields keep the usual ``dgmodel/dgimg/dgweight/dgsigma`` names;
    independent fields add ``_independent``. Both stages' coefficients and
    convergence diagnostics are retained. Sigma is residual spread, separate
    from measurement variance and fitted-background uncertainty.

    Accept a Dataset or a sensor-keyed mapping. ``pointing_correction=True``
    requires a mapping containing WIC, SI12 and SI13, loaded using
    ``read_idl(..., pointing_metadata=True)``. Shared metadata correction is
    applied once before fitting. ``pointing_options`` configures that step.
    Neither pointing correction nor the new method changes1D defaults.
    Numerical/iteration failure raises instead of returning partial fits.
    """
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
        options = dict(
            sKnots=sKnots, sine_knots=sine_knots, temporal_knot_spacing=temporal_knot_spacing, tOrder=tOrder,
            dgsigma_freeze_iteration=dgsigma_freeze_iteration, stop=stop, dampingVal=dampingVal, tukeyVal=tukeyVal,
            variance_weighted=variance_weighted, minimum_frame_quality=minimum_frame_quality, max_iterations=max_iterations
        )
        return {sensor: backgroundmodel_BS2D(data, **options) for sensor, data in sensors.items()}
    if pointing_correction:
        raise ValueError('Shared pointing correction requires a WIC/SI12/SI13 mapping')
    if isinstance(max_iterations, (bool, np.bool_)) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 2:
        raise ValueError('max_iterations must be an integer >=2')
    if isinstance(dgsigma_freeze_iteration, (bool, np.bool_)) or not isinstance(dgsigma_freeze_iteration, (int, np.integer)) or not 1 <= dgsigma_freeze_iteration < max_iterations:
        raise ValueError('dgsigma_freeze_iteration must be a positive integer smaller than max_iterations')
    if isinstance(tOrder, (bool, np.bool_)) or not isinstance(tOrder, (int, np.integer)) or tOrder < 0:
        raise ValueError('tOrder must be a nonnegative integer')
    if not isinstance(variance_weighted, (bool, np.bool_)):
        raise ValueError('variance_weighted must be boolean')
    for name, value, zero_allowed in [('stop', stop, False), ('tukeyVal', tukeyVal, False), ('dampingVal', dampingVal, True)]:
        if isinstance(value, (bool, np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value < 0 or (value == 0 and not zero_allowed):
            raise ValueError(f'{name} must be finite and {"nonnegative" if zero_allowed else "positive"}')
    if minimum_frame_quality is not None and (isinstance(minimum_frame_quality, (bool, np.bool_)) or not isinstance(minimum_frame_quality, (int, np.integer)) or minimum_frame_quality not in [0, 1, 2]):
        raise ValueError('minimum_frame_quality must be0,1,2 or None')
    required = ['img', 'sza', 'dza', 'solar_view_sine_factor'] + (['img_variance'] if variance_weighted else [])
    missing = [name for name in required if name not in imgs]
    if missing:
        raise ValueError('2D background requires fields: '+', '.join(missing))
    data = imgs.copy(deep=True)
    data = data.drop_vars([name for name in data.variables if name.startswith('dg')])
    if 'date' not in data.img.dims:
        if 'date' not in data.coords or data.date.ndim != 0:
            raise ValueError('A single image requires its scalar date coordinate')
        data = data.expand_dims('date')
    if data.img.dims != ('date', 'row', 'col') or any(data[name].dims != data.img.dims for name in required):
        raise ValueError('Image, geometry and variance fields must have dimensions(date,row,col)')
    if minimum_frame_quality is not None and 'frame_quality' in data:
        if data.frame_quality.dims != ('date',):
            raise ValueError('frame_quality must have only the date dimension')
        data = data.isel(date=data.frame_quality.values >= minimum_frame_quality)
    if not data.sizes['date']:
        raise ValueError('No frames satisfy the requested frame quality')
    if not np.issubdtype(data.date.dtype, np.datetime64) or np.any(np.isnat(data.date.values)):
        raise ValueError('Image dates must be finite')
    if sKnots is None:
        inner = np.cos(np.deg2rad([70, 80, 85, 90, 95, 100, 110, 120]))
        inner[np.abs(inner) < 1e-12] = 0
        sKnots = np.sort(np.r_[-3.5, inner, 1.5, 3.5])
    sine_knots = [-1., 1.] if sine_knots is None else sine_knots
    x_spline, sine_spline = _spline(sKnots, 'sKnots'), _spline(sine_knots, 'sine_knots')
    if sine_spline.t[0] != -1 or sine_spline.t[-1] != 1:
        raise ValueError('sine_knots must span the physical domain[-1,1]')
    minutes = (data.date.values-data.date.values[0])/np.timedelta64(1, 'm')
    time_breaks, temporal, n_time, order, spacing = bg._temporal_spline_basis(minutes, temporal_knot_spacing, tOrder)
    _, directional, _, _, _ = bg._temporal_spline_basis(minutes, temporal_knot_spacing/2, tOrder)
    ridge = _reference_ridge(data, temporal, directional, dampingVal, variance_weighted)
    temporal_fit = _fit_sequence(data, x_spline, sine_spline, temporal, ridge, dgsigma_freeze_iteration, stop, tukeyVal, max_iterations, variance_weighted)
    for field in ['model', 'img', 'weight', 'sigma']:
        data['dg'+field] = (data.img.dims, temporal_fit[field])
        data['dg'+field+'_independent'] = xr.full_like(data['dg'+field], np.nan)
    controls = ['dg_x_control', 'dg_sine_control', 'dg_time_control']
    data['dgmodel_coefficients'] = (controls, temporal_fit['coefficients'])
    data['dgsigma_fit_coefficients'] = (['dg_sigma_fit', 'dg_sigma_parameter'], temporal_fit['sigma_history'])
    data = data.assign_coords(dg_sigma_fit=np.arange(1, dgsigma_freeze_iteration+1), dg_sigma_parameter=['a', 'b'])
    data['dgfit_iterations'] = temporal_fit['iterations']
    data['dgfit_converged'] = True
    data['dgfit_change'] = temporal_fit['change']
    data['dgsolver_flag'] = 1
    sigma_law = temporal_fit['sigma_coefficients']

    # Full frame backgrounds; never inherit the temporal shape or pixel mask.
    coefficients, iterations, changes, ridges = [], [], [], []
    for frame in range(data.sizes['date']):
        image = data.isel(date=[frame])
        singleton = np.ones((1, 1))
        local_ridge = _reference_ridge(image, singleton, singleton, dampingVal, variance_weighted)
        fitted = _fit_sequence(image, x_spline, sine_spline, singleton, local_ridge, dgsigma_freeze_iteration, stop, tukeyVal, max_iterations, variance_weighted, sigma_coefficients=sigma_law)
        for field in ['model', 'img', 'weight', 'sigma']:
            data['dg'+field+'_independent'].values[frame] = fitted[field][0]
        coefficients.append(fitted['coefficients'][:, :, 0])
        iterations.append(fitted['iterations'])
        changes.append(fitted['change'])
        ridges.append(local_ridge)
    data['dgmodel_independent_coefficients'] = (['date', *controls[:2]], np.asarray(coefficients))
    data['dgfit_iterations_independent'] = ('date', iterations)
    data['dgfit_converged_independent'] = ('date', np.ones(data.sizes['date'], dtype=bool))
    data['dgfit_change_independent'] = ('date', changes)
    data['dg_ridge_independent'] = ('date', ridges)
    data['dgsolver_flag_independent'] = ('date', np.ones(data.sizes['date'], dtype=np.int8))
    common = {
        'spatial_model': 'x_sine_2d', 'spatial_knots': np.asarray(sKnots), 'sine_knots': np.asarray(sine_knots),
        'spatial_order': 3, 'spatial_basis_scale': float(np.sqrt(sine_spline.c.shape[0])),
        'sigma_model_coefficients': sigma_law, 'dgsigma_freeze_iteration': int(dgsigma_freeze_iteration),
        'monotonic_constraint': 'nondecreasing_x_at_every_sine_time_control', 'tukey_value': float(tukeyVal),
        'convergence_tolerance': float(stop), 'measurement_variance': 'img_variance' if variance_weighted else 'not_used'
    }
    data.dgmodel.attrs = {
        **common, 'long_name': 'Temporal2D background', 'noise_model': 'binned',
        'dgsigma_estimation': 'pooled_per_fit_then_frozen_pixel_field',
        'temporal_knots_minutes_from_first_image': time_breaks, 'temporal_order': int(order),
        'temporal_knot_spacing_target_minutes': float(temporal_knot_spacing),
        'temporal_knot_spacing_actual_minutes': float(spacing), 'ridge_scalar': ridge,
        'ridge_reference': 'legacy coarse base plus tapered sine; reference design only sets scalar damping'
    }
    data.dgmodel_independent.attrs = {
        **common, 'long_name': 'Independent2D background', 'temporal_order': 0,
        'noise_model': 'supplied_quadratic', 'dgsigma_estimation': 'supplied_relation_each_iteration',
        'robust_initialization': 'fresh all-ones; updated every fit', 'shape_prior': 'none'
    }
    for suffix in ['', '_independent']:
        data['dgsigma'+suffix].attrs = {'long_name': 'Empirical residual spread; not detector noise or background-fit uncertainty'}
        data['dgimg'+suffix].attrs = {'long_name': 'Counts minus '+('independent' if suffix else 'temporal')+'2D background'}
        data['dgweight'+suffix].attrs = {'long_name': 'Final Tukey bisquare objective weights'}
    data.dgmodel_coefficients.attrs = {'ordering': 'x,sine,time', 'basis_scale': common['spatial_basis_scale']}
    data.dgmodel_independent_coefficients.attrs = {'ordering': 'date,x,sine', 'basis_scale': common['spatial_basis_scale']}
    return data
