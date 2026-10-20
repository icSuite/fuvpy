"""Joint metadata-only IMAGE/FUV pointing repair on each camera's native grid.

The WIC time series supplies one spacecraft correction shared by all cameras.
Camera alignments, brightness, and detector variance never determine a repair.
The historical emission ellipsoid, SUN3 and 1997 Apex grid are retained.
"""

from collections.abc import Mapping
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates
from scipy.spatial.transform import Rotation

from .viewing_geometry import calculate_viewing_geometry
from ..utils.sunlight import subsol


SENSORS = {'WIC', 'SI12', 'SI13'}
GEOGRAPHIC_FIELDS = ['glat', 'glon', 'sza', 'dza']
MAGNETIC_FIELDS = ['mlat', 'mlon', 'mlt']
SOLAR_FIELDS = ['solar_view_cosine_factor', 'solar_view_sine_factor']

ROBUST_STATE_DEFAULTS = {
    'maximum_gap_seconds': 300.,
    'simultaneous_tolerance_seconds': 1.,
    'minimum_span_observations': 5,
    'prediction_window': 11,
    'outlier_sigma': 8.,
    'position_outlier_floor_km': 20.,
    'axis_outlier_floor_deg': .05,
    'phase_outlier_floor_deg': .25,
    'position_measurement_sigma_km': 1.,
    'axis_measurement_sigma_deg': .01,
    'phase_measurement_sigma_deg': .02,
    'position_acceleration_sigma_km_s2': .01,
    'axis_acceleration_sigma_deg_s2': 1e-5,
    'phase_acceleration_sigma_deg_s2': 1e-5,
    'maximum_position_uncertainty_km': 25.,
    'maximum_axis_uncertainty_deg': .2,
    'maximum_phase_uncertainty_deg': 1.,
    'maximum_coordinate_displacement_deg': 15.
}


def _historical_sun(time):
    """FUVVIEW3 SUN3 inertial Sun, including its single-precision literals."""
    time = pd.Timestamp(time)
    value = lambda x: float(np.float32(x))
    pi = value(np.pi)
    milliseconds = (time.hour*3600 + time.minute*60 + time.second)*1000 + time.microsecond//1000
    seconds = np.float32(milliseconds)/np.float32(1000)
    dj = 365*(time.year - 1900) + (time.year - 1901)//4 + time.dayofyear + float(seconds/np.float32(86400)) - .5
    t = dj/36525
    vl = (value(279.696678) + value(.9856473354)*dj) % 360
    g = (value(358.475845) + value(.985600267)*dj) % 360
    slon = vl + (value(1.91946) - value(.004789)*t)*np.sin(g*pi/180) + value(.020094)*np.sin(2*g*pi/180)
    obliquity = value(23.45229) - value(.0130125)*t
    slp = slon - value(.005686)
    sine = np.sin(obliquity*pi/180)*np.sin(slp*pi/180)
    cosine = np.sqrt(1 - sine*sine)
    dec = np.arctan(sine/cosine)*180/pi
    cotob = np.cos(obliquity*pi/180)/np.sin(obliquity*pi/180)
    ra = 180 - np.arctan2(cotob*sine/cosine, -np.cos(slp*pi/180)/cosine)*180/pi
    ra, dec = np.deg2rad([ra, dec])
    return np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])


def _ellipsoid(glat, glon, height):
    latitude, longitude = np.deg2rad([glat, glon])
    equator = 6378.136 + height
    pole = 6378.136*(1 - 1/298.257) + height
    ratio = (pole/equator)**2
    radius = equator/np.sqrt(1 - (1 - ratio)*np.sin(latitude)**2)
    point = np.stack([radius*np.cos(latitude)*np.cos(longitude), radius*np.cos(latitude)*np.sin(longitude), radius*ratio*np.sin(latitude)], axis=-1)
    return point, equator, ratio


def _triad(observer, sun):
    first = observer/np.linalg.norm(observer)
    third = np.cross(first, sun)
    norm = np.linalg.norm(third)
    if norm < 1e-3:
        raise ValueError('Observer and Sun do not constrain the Earth-frame rotation')
    third /= norm
    return np.column_stack([first, np.cross(third, first), third])


def _native_geometry(glat, glon, sza, dza, observer, centre_latitude, centre_longitude, height, time):
    """Recover rays from coordinates; this also respects downshifted SI grids.

    GCI-to-geographic rotation includes historical precession/nutation, so a
    z rotation inferred from observer longitude alone is insufficient. Stored
    SZA constrains the geographic Sun; the observer supplies a second vector.
    """
    finite = np.isfinite(glat) & np.isfinite(glon) & np.isfinite(sza) & np.isfinite(dza)
    valid = finite & (np.abs(glat) <= 90) & (dza >= 0) & (dza < 75) & (sza >= 0) & (sza <= 180)
    if valid.sum() < 6:
        raise ValueError('At least six valid native coordinate pixels are needed to reconstruct pointing')
    safe_latitude = np.where(finite & (np.abs(glat) <= 90), glat, np.nan)
    safe_longitude = np.where(finite & (np.abs(glat) <= 90), glon, np.nan)
    point, equator, ratio = _ellipsoid(safe_latitude, safe_longitude, height)
    up = point/np.linalg.norm(point, axis=-1, keepdims=True)
    if np.linalg.cond(up[valid]) > 1e6:
        raise ValueError('Native coordinate pixels do not constrain the Earth-frame rotation')
    sun = _historical_sun(time)
    geographic_sun = np.linalg.lstsq(up[valid], np.cos(np.deg2rad(sza[valid])), rcond=None)[0]
    geographic_sun /= np.linalg.norm(geographic_sun)
    latitude, longitude = np.deg2rad([centre_latitude, centre_longitude])
    geographic_observer = np.linalg.norm(observer)*np.array([np.cos(latitude)*np.cos(longitude), np.cos(latitude)*np.sin(longitude), np.sin(latitude)])
    transform = _triad(geographic_observer, geographic_sun) @ _triad(observer, sun).T
    inertial = point @ transform
    look = inertial - observer
    look /= np.linalg.norm(look, axis=-1, keepdims=True)
    radial = inertial/np.linalg.norm(inertial, axis=-1, keepdims=True)
    reconstructed_sza = np.rad2deg(np.arccos(np.clip(radial @ sun, -1, 1)))
    reconstructed_dza = np.rad2deg(np.arccos(np.clip(np.sum(radial*(-look), axis=-1), -1, 1)))
    error = max(np.max(np.abs(reconstructed_sza[valid] - sza[valid])), np.max(np.abs(reconstructed_dza[valid] - dza[valid])))
    if not np.isfinite(error) or error > 1e-3:
        raise ValueError(f'Native pointing reconstruction disagrees with stored angles by {error:.6g} degrees')
    return look, transform, sun, equator, ratio, float(error)


@lru_cache(maxsize=1)
def _apex_grid():
    filename = Path(__file__).resolve().parents[1]/'data/pointing/mlatlon.1997a.xdr'
    if hashlib.sha256(filename.read_bytes()).hexdigest() != '54894228da8b60dc21cc5337cf112c4f9e40c26c65d308cf31034562f4e5b549':
        raise ValueError('Historical Apex lookup grid does not match its recorded source checksum')
    values = np.fromfile(filename, dtype='>f4').astype(float)
    if values.size != 2*181*361:
        raise ValueError('Invalid historical Apex lookup grid')
    latitude = values[:181*361].reshape(181, 361)
    longitude = np.deg2rad(values[181*361:].reshape(181, 361))
    return latitude, np.sin(longitude), np.cos(longitude)


def _magnetic(glat, glon):
    """Same 1997 lookup interpolation as FUVVIEW3 geo_to_apex.pro."""
    latitude, sine, cosine = _apex_grid()
    coordinates = [np.asarray(glat) + 90, np.mod(glon, 360)]
    mlat = map_coordinates(latitude, coordinates, order=1, mode='nearest')
    mlon = np.rad2deg(np.arctan2(map_coordinates(sine, coordinates, order=1, mode='nearest'), map_coordinates(cosine, coordinates, order=1, mode='nearest')))
    return mlat, mlon


def _local_predictions(time, values, segment, window, degree, omit_centre=False, robust=False):
    prediction = values.copy()
    supported = np.zeros(len(time), dtype=bool)
    polynomials = [None]*len(time)
    for group in np.unique(segment):
        indices = np.flatnonzero(segment == group)
        if len(indices) < window:
            continue
        for offset, index in enumerate(indices):
            lower = max(0, min(offset - window//2, len(indices) - window))
            neighbours = indices[lower:lower + window]
            if omit_centre:
                neighbours = neighbours[neighbours != index]
            x = (time[neighbours] - time[index])/600
            design = np.vander(x, degree + 1, increasing=True)
            delta = values[neighbours] - values[index]
            weights = np.ones(len(neighbours))
            for _ in range(8 if robust else 1):
                coefficients = np.linalg.lstsq(design*weights[:, None], delta*weights[:, None], rcond=None)[0]
                residual = np.linalg.norm(delta - design @ coefficients, axis=1)
                scale = max(1e-6, 1.4826*np.median(np.abs(residual - np.median(residual))))
                weights = np.minimum(1, 1.5*scale/np.maximum(residual, 1e-12))
            prediction[index] = values[index] + coefficients[0]
            supported[index] = True
            polynomials[index] = coefficients
    return prediction, supported, polynomials


def _reference_corrections(images, window, position_floor):
    time = (images.pointing_time.values.astype('datetime64[ns]').astype(np.int64) - images.pointing_time.values[0].astype('datetime64[ns]').astype(np.int64))/1e9
    spin = images.SPIN.values.astype(np.int64)
    intervals = np.diff(time)
    cadence = np.median(intervals) if len(intervals) else np.inf
    segment = np.cumsum(np.r_[True, (intervals > 1.5*cadence) | (np.diff(spin) != 1)])
    position = np.asarray(images.O_GCI.values, dtype=float)
    prediction, supported, polynomials = _local_predictions(time, position, segment, window, 3, True, True)
    residual = np.linalg.norm(position - prediction, axis=1)
    centre = np.median(residual[supported]) if np.any(supported) else 0.
    scale = 1.4826*np.median(np.abs(residual[supported] - centre)) if np.any(supported) else 0.
    threshold = max(position_floor, centre + 10*scale)
    outlier = np.zeros(len(time), dtype=bool)
    for index in range(2, len(time) - 2):
        neighbours = residual[index - 2:index + 3].copy()
        neighbours[2] = 0
        outlier[index] = supported[index] and residual[index] > threshold and residual[index] > 1.5*neighbours.max()
    phase = np.rad2deg(np.unwrap(np.deg2rad(images.SPIN_PHASE.values))) + 360*(spin - spin[0])
    filtered, _, _ = _local_predictions(time, phase[:, None], segment, window, 2)
    return {'time': time, 'spin': spin, 'segment': segment, 'position': position, 'outlier': outlier, 'polynomials': polynomials, 'phase_delta': filtered[:, 0] - phase, 'supported': supported, 'threshold': threshold, 'residual': residual}


def _state_settings(options):
    settings = dict(ROBUST_STATE_DEFAULTS)
    if options is not None:
        if not isinstance(options, dict):
            raise TypeError('state_options must be a dictionary or None')
        unknown = set(options) - set(settings)
        if unknown:
            raise ValueError('Unknown robust-state option(s): ' + ', '.join(sorted(unknown)))
        settings.update(options)

    integers = ['minimum_span_observations', 'prediction_window']
    for name in integers:
        if isinstance(settings[name], (bool, np.bool_)) or not isinstance(settings[name], (int, np.integer)):
            raise ValueError(f'{name} must be an integer')
    if settings['minimum_span_observations'] < 3:
        raise ValueError('minimum_span_observations must be at least three')
    if settings['prediction_window'] < 5 or settings['prediction_window'] % 2 != 1:
        raise ValueError('prediction_window must be an odd integer of at least five')

    for name, value in settings.items():
        if name in integers:
            continue
        if isinstance(value, (bool, np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be finite and positive')
    return settings


def _state_segments(time, spin, maximum_gap_seconds):
    """Split WIC records at time gaps and counter discontinuities."""
    if len(time) == 0:
        return np.empty(0, dtype=int)
    breaks = np.r_[True, (np.diff(time) > maximum_gap_seconds) | (np.diff(spin) != 1)]
    return np.cumsum(breaks) - 1


def _unwrapped_phase(phase, spin, segment):
    result = np.full(len(phase), np.nan)
    for group in np.unique(segment):
        indices = np.flatnonzero(segment == group)
        local = np.rad2deg(np.unwrap(np.deg2rad(phase[indices])))
        result[indices] = local + 360*(spin[indices] - spin[indices][0])
    return result


def _prediction_threshold(residual, supported, segment, floor, sigma):
    threshold = np.full(len(residual), np.inf)
    for group in np.unique(segment):
        indices = np.flatnonzero((segment == group) & supported & np.isfinite(residual))
        if len(indices) == 0:
            continue
        centre = np.median(residual[indices])
        scale = 1.4826*np.median(np.abs(residual[indices] - centre))
        threshold[segment == group] = max(floor, centre + sigma*scale)
    return threshold


def _isolated_candidates(candidate, segment):
    """Require two ordinary observations on both sides of an outlier."""
    isolated = np.zeros(len(candidate), dtype=bool)
    for group in np.unique(segment):
        indices = np.flatnonzero(segment == group)
        for offset in range(2, len(indices) - 2):
            index = indices[offset]
            neighbours = indices[offset - 2:offset + 3]
            isolated[index] = candidate[index] and not np.any(candidate[neighbours[[0, 1, 3, 4]]])
    return isolated


def _split_around_unresolved(segment, unresolved):
    """Prevent a smoother from joining observations across ambiguous records."""
    result = np.full(len(segment), -1, dtype=int)
    group_number = 0
    for group in np.unique(segment):
        indices = np.flatnonzero(segment == group)
        run = []
        for index in indices:
            if unresolved[index]:
                if run:
                    result[run] = group_number
                    group_number += 1
                    run = []
                continue
            run.append(index)
        if run:
            result[run] = group_number
            group_number += 1
    return result


def _constant_rate_smoother(time, values, measurement_sigma, acceleration_sigma, weights):
    """Linear constant-rate Kalman filter followed by an RTS smoother."""
    values = np.asarray(values, dtype=float)
    dimensions = values.shape[1]
    count = len(time)
    identity = np.eye(dimensions)
    observation = np.column_stack([identity, np.zeros_like(identity)])
    median_step = np.median(np.diff(time)) if count > 1 else 1.
    initial_rate = (values[-1] - values[0])/max(time[-1] - time[0], median_step)
    state = np.r_[values[0], initial_rate]
    covariance = np.diag(np.r_[np.full(dimensions, (10*measurement_sigma)**2), np.full(dimensions, (10*measurement_sigma/median_step)**2)])

    filtered_state = np.empty((count, 2*dimensions))
    filtered_covariance = np.empty((count, 2*dimensions, 2*dimensions))
    predicted_state = np.empty_like(filtered_state)
    predicted_covariance = np.empty_like(filtered_covariance)
    transitions = [None]*count

    for index in range(count):
        if index:
            step = time[index] - time[index - 1]
            transition = np.block([[identity, step*identity], [np.zeros_like(identity), identity]])
            process = acceleration_sigma**2*np.block([[(step**4/4)*identity, (step**3/2)*identity], [(step**3/2)*identity, step**2*identity]])
            state = transition @ state
            covariance = transition @ covariance @ transition.T + process
            transitions[index] = transition
        predicted_state[index] = state
        predicted_covariance[index] = covariance

        noise = (measurement_sigma/max(float(weights[index]), 1e-4))**2*identity
        innovation_covariance = observation @ covariance @ observation.T + noise
        gain = np.linalg.solve(innovation_covariance, observation @ covariance).T
        state = state + gain @ (values[index] - observation @ state)
        covariance = (np.eye(2*dimensions) - gain @ observation) @ covariance
        covariance = (covariance + covariance.T)/2
        filtered_state[index] = state
        filtered_covariance[index] = covariance

    smoothed_state = filtered_state.copy()
    smoothed_covariance = filtered_covariance.copy()
    for index in range(count - 2, -1, -1):
        transition = transitions[index + 1]
        gain = np.linalg.solve(predicted_covariance[index + 1], transition @ filtered_covariance[index]).T
        smoothed_state[index] += gain @ (smoothed_state[index + 1] - predicted_state[index + 1])
        smoothed_covariance[index] += gain @ (smoothed_covariance[index + 1] - predicted_covariance[index + 1]) @ gain.T
        smoothed_covariance[index] = (smoothed_covariance[index] + smoothed_covariance[index].T)/2
    return smoothed_state, smoothed_covariance


def _robust_wic_state(images, settings):
    """Estimate position, spin axis and phase from WIC metadata only."""
    epoch = images.pointing_time.values[0].astype('datetime64[ns]').astype(np.int64)
    time = (images.pointing_time.values.astype('datetime64[ns]').astype(np.int64) - epoch)/1e9
    spin = images.SPIN.values.astype(np.int64)
    position = np.asarray(images.O_GCI.values, dtype=float)
    axis = np.asarray(images.A_GCI3.values, dtype=float)
    axis /= np.linalg.norm(axis, axis=1, keepdims=True)
    segment = _state_segments(time, spin, settings['maximum_gap_seconds'])
    phase = _unwrapped_phase(np.asarray(images.SPIN_PHASE.values, dtype=float), spin, segment)

    position_prediction, position_supported, _ = _local_predictions(time, position, segment, settings['prediction_window'], 3, True, True)
    axis_prediction, axis_supported, _ = _local_predictions(time, axis, segment, settings['prediction_window'], 2, True, True)
    phase_prediction, phase_supported, _ = _local_predictions(time, phase[:, None], segment, settings['prediction_window'], 2, True, True)
    axis_prediction /= np.linalg.norm(axis_prediction, axis=1, keepdims=True)

    position_residual = np.linalg.norm(position - position_prediction, axis=1)
    axis_residual = np.rad2deg(np.arccos(np.clip(np.sum(axis*axis_prediction, axis=1), -1, 1)))
    phase_residual = np.abs(phase - phase_prediction[:, 0])
    position_threshold = _prediction_threshold(position_residual, position_supported, segment, settings['position_outlier_floor_km'], settings['outlier_sigma'])
    axis_threshold = _prediction_threshold(axis_residual, axis_supported, segment, settings['axis_outlier_floor_deg'], settings['outlier_sigma'])
    phase_threshold = _prediction_threshold(phase_residual, phase_supported, segment, settings['phase_outlier_floor_deg'], settings['outlier_sigma'])
    candidate = ((position_supported & (position_residual > position_threshold)) |
                 (axis_supported & (axis_residual > axis_threshold)) |
                 (phase_supported & (phase_residual > phase_threshold)))
    isolated = _isolated_candidates(candidate, segment)
    unresolved = candidate & ~isolated
    fit_segment = _split_around_unresolved(segment, unresolved)

    fitted_position = np.full_like(position, np.nan)
    fitted_velocity = np.full_like(position, np.nan)
    fitted_axis = np.full_like(axis, np.nan)
    fitted_axis_rate = np.full_like(axis, np.nan)
    fitted_phase = np.full_like(phase, np.nan)
    fitted_phase_rate = np.full_like(phase, np.nan)
    position_uncertainty = np.full(len(time), np.inf)
    axis_uncertainty = np.full(len(time), np.inf)
    phase_uncertainty = np.full(len(time), np.inf)
    supported = np.zeros(len(time), dtype=bool)
    weights = np.where(isolated, .01, 1.)

    for group in np.unique(fit_segment[fit_segment >= 0]):
        indices = np.flatnonzero(fit_segment == group)
        if len(indices) < settings['minimum_span_observations']:
            unresolved[indices] = True
            continue
        local_time = time[indices] - time[indices][0]
        position_state, position_covariance = _constant_rate_smoother(local_time, position[indices], settings['position_measurement_sigma_km'], settings['position_acceleration_sigma_km_s2'], weights[indices])
        axis_state, axis_covariance = _constant_rate_smoother(local_time, axis[indices], np.deg2rad(settings['axis_measurement_sigma_deg']), np.deg2rad(settings['axis_acceleration_sigma_deg_s2']), weights[indices])
        phase_state, phase_covariance = _constant_rate_smoother(local_time, phase[indices, None], settings['phase_measurement_sigma_deg'], settings['phase_acceleration_sigma_deg_s2'], weights[indices])

        fitted_position[indices] = position_state[:, :3]
        fitted_velocity[indices] = position_state[:, 3:]
        fitted_axis[indices] = axis_state[:, :3]/np.linalg.norm(axis_state[:, :3], axis=1, keepdims=True)
        fitted_axis_rate[indices] = axis_state[:, 3:]
        fitted_phase[indices] = phase_state[:, 0]
        fitted_phase_rate[indices] = phase_state[:, 1]
        position_uncertainty[indices] = np.sqrt(np.maximum(0, np.trace(position_covariance[:, :3, :3], axis1=1, axis2=2)))
        axis_uncertainty[indices] = np.rad2deg(np.sqrt(np.maximum(0, np.trace(axis_covariance[:, :3, :3], axis1=1, axis2=2))))
        phase_uncertainty[indices] = np.sqrt(np.maximum(0, phase_covariance[:, 0, 0]))
        supported[indices] = ((position_uncertainty[indices] <= settings['maximum_position_uncertainty_km']) &
                              (axis_uncertainty[indices] <= settings['maximum_axis_uncertainty_deg']) &
                              (phase_uncertainty[indices] <= settings['maximum_phase_uncertainty_deg']))

    return {
        'epoch': epoch, 'time': time, 'segment': fit_segment, 'position': position, 'axis': axis, 'phase': phase,
        'fitted_position': fitted_position, 'fitted_velocity': fitted_velocity, 'fitted_axis': fitted_axis,
        'fitted_axis_rate': fitted_axis_rate, 'fitted_phase': fitted_phase, 'fitted_phase_rate': fitted_phase_rate,
        'position_uncertainty': position_uncertainty, 'axis_uncertainty': axis_uncertainty,
        'phase_uncertainty': phase_uncertainty, 'supported': supported, 'isolated_outlier': isolated,
        'unresolved': unresolved, 'position_threshold': position_threshold,
        'axis_threshold': axis_threshold, 'phase_threshold': phase_threshold
    }


def _validate_images(images, allow_empty=False):
    required = ['pointing_time', 'O_GCI', 'A_GCI3', 'SPIN_PHASE', 'SPIN', 'LAT_CTR', 'LON_CTR', 'EMIS_HGT'] + GEOGRAPHIC_FIELDS + MAGNETIC_FIELDS + SOLAR_FIELDS
    valid_axes = {}
    for sensor, data in images.items():
        if 'id' in data and str(data.id.item()) != sensor:
            raise ValueError(f'{sensor} mapping entry contains sensor id {data.id.item()}')
        missing = [name for name in required if name not in data]
        if missing:
            raise ValueError(f'{sensor} pointing metadata missing: {", ".join(missing)}; load with pointing_metadata=True')
        if data.sizes.get('date', 0) < 1 and not allow_empty:
            raise ValueError(f'{sensor} requires at least one frame')
        metadata_valid = np.asarray(data.pointing_valid.values, dtype=bool) if 'pointing_valid' in data else np.ones(data.sizes['date'], dtype=bool)
        time = data.pointing_time.values.astype('datetime64[ns]')
        if np.any(np.isnat(time[metadata_valid])) or np.any(np.diff(time[metadata_valid].astype(np.int64)) <= 0):
            raise ValueError(f'{sensor} pointing_time must be finite and strictly increasing')
        for name in required[1:8]:
            if not np.all(np.isfinite(data[name].values[metadata_valid])):
                raise ValueError(f'{sensor} {name} metadata must be finite')
        spin = data.SPIN.values
        if np.any(spin[metadata_valid] != np.rint(spin[metadata_valid])):
            raise ValueError(f'{sensor} SPIN must contain integers')
        norm = np.linalg.norm(data.A_GCI3.values, axis=1)
        valid_axes[sensor] = metadata_valid & (np.abs(norm - 1) <= 1e-3)
        raw_dza = data.get('pointing_DZA', data.dza).values
        viewing_limit = float(data.attrs.get('pointing_dzalim', 80))
        visible = np.isfinite(raw_dza) & (raw_dza >= 0) & (raw_dza < viewing_limit)
        has_rays = np.any(visible, axis=tuple(range(1, visible.ndim)))
        if np.any(metadata_valid & ~valid_axes[sensor] & has_rays):
            raise ValueError(f'{sensor} A_GCI3 must contain unit spin-axis vectors in frames with visible rays')
        if data.attrs.get('pointing_correction_applied', False):
            raise ValueError(f'{sensor} pointing has already been corrected')
    if not allow_empty and not np.any(valid_axes['WIC']):
        raise ValueError('WIC has no valid spin-axis frames for shared pointing')
    return valid_axes


def _valid_wic_reference(data):
    """Keep spacecraft states suitable for fitting, regardless of image quality."""
    required = ['pointing_time', 'O_GCI', 'A_GCI3', 'SPIN_PHASE', 'SPIN']
    missing = [name for name in required if name not in data]
    if missing:
        raise ValueError(f'WIC reference metadata missing: {", ".join(missing)}')
    if 'id' in data and str(data.id.item()) != 'WIC':
        raise ValueError('Shared pointing reference must be WIC')
    if data.attrs.get('pointing_correction_applied', False):
        raise ValueError('WIC pointing reference has already been corrected')

    time = data.pointing_time.values.astype('datetime64[ns]')
    position = np.asarray(data.O_GCI.values, dtype=float)
    axis = np.asarray(data.A_GCI3.values, dtype=float)
    phase = np.asarray(data.SPIN_PHASE.values, dtype=float)
    spin = np.asarray(data.SPIN.values, dtype=float)
    valid = ~np.isnat(time) & np.all(np.isfinite(position), axis=1) & np.all(np.isfinite(axis), axis=1)
    valid &= np.isfinite(phase) & np.isfinite(spin) & (spin == np.rint(spin))
    valid &= np.abs(np.linalg.norm(axis, axis=1) - 1) <= 1e-3
    if 'pointing_valid' in data:
        valid &= np.asarray(data.pointing_valid.values, dtype=bool)
    if not np.any(valid):
        raise ValueError('WIC has no valid spacecraft states for shared pointing')
    if np.any(np.diff(time[valid].astype(np.int64)) <= 0):
        raise ValueError('WIC reference pointing_time must be strictly increasing')
    return valid


def _correct_frame(data, index, observer, phase_delta=0., correction_rotation=None):
    raw = {name: np.asarray(data.get('pointing_' + name.upper(), data[name]).isel(date=index).values, dtype=float) for name in GEOGRAPHIC_FIELDS}
    height = float(data.EMIS_HGT.values[index])
    time = data.pointing_time.values[index]
    original_observer = np.asarray(data.O_GCI.values[index], dtype=float)
    look, transform, sun, equator, ratio, error = _native_geometry(raw['glat'], raw['glon'], raw['sza'], raw['dza'], original_observer, float(data.LAT_CTR.values[index]), float(data.LON_CTR.values[index]), height, time)
    if correction_rotation is None:
        axis = np.asarray(data.A_GCI3.values[index], dtype=float).copy()
        axis /= np.linalg.norm(axis)
        correction_rotation = Rotation.from_rotvec(np.deg2rad(phase_delta)*axis).as_matrix()
    look = look @ correction_rotation.T
    metric = np.array([1, 1, 1/ratio])
    aa = np.sum(metric*look**2, axis=-1)
    bb = 2*np.sum(metric*look*observer, axis=-1)
    discriminant = bb*bb - 4*aa*(np.sum(metric*observer**2) - equator**2)
    distance = (-bb - np.sqrt(np.maximum(discriminant, 0)))/(2*aa)
    point = observer + distance[..., None]*look
    radial = point/np.linalg.norm(point, axis=-1, keepdims=True)
    dza = np.rad2deg(np.arccos(np.clip(np.sum(radial*(-look), axis=-1), -1, 1)))
    sza = np.rad2deg(np.arccos(np.clip(radial @ sun, -1, 1)))
    geographic = point @ transform.T
    glat = np.rad2deg(np.arctan2(geographic[..., 2], ratio*np.hypot(geographic[..., 0], geographic[..., 1])))
    glon = np.rad2deg(np.arctan2(geographic[..., 1], geographic[..., 0])) % 360
    geographic_observer = observer @ transform.T
    centre_latitude = np.rad2deg(np.arctan2(geographic_observer[2], np.hypot(geographic_observer[0], geographic_observer[1])))
    centre_longitude = np.rad2deg(np.arctan2(geographic_observer[1], geographic_observer[0])) % 360
    subsolar_latitude, subsolar_longitude = subsol([pd.Timestamp(data.date.values[index])])
    view = calculate_viewing_geometry(glat, glon, np.linalg.norm(observer), centre_latitude, centre_longitude, height, subsolar_latitude[0], subsolar_longitude[0])
    support = (discriminant > 0) & (distance > 0) & (dza < float(data.attrs.get('pointing_dzalim', 80)))
    for name in GEOGRAPHIC_FIELDS + MAGNETIC_FIELDS:
        support &= np.isfinite(data[name].values[index])
    if np.any(support) and np.max(np.abs(view['dza'][support] - dza[support])) > 1e-3:
        raise ValueError('Corrected ray and geographic viewing directions disagree')
    baseline = np.isfinite(data.mlat.values[index]) & np.isfinite(data.mlon.values[index]) & np.isfinite(data.mlt.values[index])
    offset = (data.mlon.values[index][baseline] - 15*data.mlt.values[index][baseline] + 180) % 360 - 180
    centre_offset = np.rad2deg(np.angle(np.mean(np.exp(1j*np.deg2rad(offset))))) if len(offset) else np.nan
    if len(offset) and np.max(np.abs((offset - centre_offset + 180) % 360 - 180)) > 1e-3:
        raise ValueError('Stored magnetic longitude and local time do not share one subsolar offset')
    mlat, mlon = _magnetic(glat, glon)
    fields = {'glat': glat, 'glon': glon, 'sza': sza, 'dza': dza, 'mlat': mlat, 'mlon': mlon, 'mlt': ((mlon - centre_offset)/15) % 24}
    fields.update({name: view[name] for name in SOLAR_FIELDS})
    fields = {name: np.where(support, values, np.nan) for name, values in fields.items()}
    return fields, centre_latitude, centre_longitude, error


def _attitude(axis, phase, reference=None):
    """Construct a deterministic inertial attitude from spin axis and phase."""
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    if reference is None:
        reference = np.eye(3)[np.argmin(np.abs(axis))]
    reference = np.asarray(reference, dtype=float)
    first = reference - np.dot(reference, axis)*axis
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    zero_phase = np.column_stack([first, second, axis])
    spin = Rotation.from_rotvec(np.deg2rad(np.mod(phase, 360))*axis).as_matrix()
    return spin @ zero_phase


def _interpolate_wic_state(reference, target_time, tolerance):
    """Interpolate one supported WIC state without crossing a span boundary."""
    for group in np.unique(reference['segment'][reference['supported']]):
        supported_indices = np.flatnonzero((reference['segment'] == group) & reference['supported'])
        runs = np.split(supported_indices, np.flatnonzero(np.diff(supported_indices) != 1) + 1)
        for indices in runs:
            if len(indices) == 0:
                continue
            first_time, last_time = reference['time'][indices[[0, -1]]]
            if target_time < first_time - tolerance or target_time > last_time + tolerance:
                continue

            clipped_time = np.clip(target_time, first_time, last_time)
            interpolate = lambda values: np.array([np.interp(clipped_time, reference['time'][indices], values[indices, component]) for component in range(values.shape[1])])
            raw_position = interpolate(reference['position'])
            raw_axis = interpolate(reference['axis'])
            raw_axis /= np.linalg.norm(raw_axis)
            raw_phase = np.interp(clipped_time, reference['time'][indices], reference['phase'][indices])
            position = interpolate(reference['fitted_position'])
            axis = interpolate(reference['fitted_axis'])
            axis /= np.linalg.norm(axis)
            phase = np.interp(clipped_time, reference['time'][indices], reference['fitted_phase'][indices])

            if target_time != clipped_time:
                endpoint = indices[0] if target_time < first_time else indices[-1]
                step = target_time - reference['time'][endpoint]
                position = reference['fitted_position'][endpoint] + step*reference['fitted_velocity'][endpoint]
                axis = reference['fitted_axis'][endpoint] + step*reference['fitted_axis_rate'][endpoint]
                axis /= np.linalg.norm(axis)
                phase = reference['fitted_phase'][endpoint] + step*reference['fitted_phase_rate'][endpoint]
                if len(indices) > 1:
                    neighbour = indices[1] if endpoint == indices[0] else indices[-2]
                    raw_step = reference['time'][endpoint] - reference['time'][neighbour]
                    raw_position = reference['position'][endpoint] + step*(reference['position'][endpoint] - reference['position'][neighbour])/raw_step
                    raw_axis = reference['axis'][endpoint] + step*(reference['axis'][endpoint] - reference['axis'][neighbour])/raw_step
                    raw_axis /= np.linalg.norm(raw_axis)
                    raw_phase = reference['phase'][endpoint] + step*(reference['phase'][endpoint] - reference['phase'][neighbour])/raw_step

            position_uncertainty = np.interp(clipped_time, reference['time'][indices], reference['position_uncertainty'][indices])
            axis_uncertainty = np.interp(clipped_time, reference['time'][indices], reference['axis_uncertainty'][indices])
            phase_uncertainty = np.interp(clipped_time, reference['time'][indices], reference['phase_uncertainty'][indices])
            nearest = indices[np.argmin(np.abs(reference['time'][indices] - clipped_time))]
            return {
                'raw_position': raw_position, 'raw_axis': raw_axis, 'raw_phase': raw_phase,
                'position': position, 'axis': axis, 'phase': phase,
                'position_uncertainty': position_uncertainty, 'axis_uncertainty': axis_uncertainty,
                'phase_uncertainty': phase_uncertainty, 'isolated_outlier': reference['isolated_outlier'][nearest]
            }
    return None


def _maximum_geographic_displacement(original, corrected):
    latitude1 = np.deg2rad(original['glat'])
    latitude2 = np.deg2rad(corrected['glat'])
    longitude1 = np.deg2rad(original['glon'])
    longitude2 = np.deg2rad(corrected['glon'])
    common = np.isfinite(latitude1) & np.isfinite(latitude2) & np.isfinite(longitude1) & np.isfinite(longitude2)
    if not np.any(common):
        return np.inf
    cosine = np.sin(latitude1[common])*np.sin(latitude2[common]) + np.cos(latitude1[common])*np.cos(latitude2[common])*np.cos(longitude2[common] - longitude1[common])
    return float(np.max(np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))))


def _correct_pointing_robust_state(images, reference_wic, settings):
    valid_axes = _validate_images(images, allow_empty=reference_wic is not None)
    if reference_wic is None:
        reference_wic = images['WIC']
        reference_valid = valid_axes['WIC']
    else:
        reference_valid = _valid_wic_reference(reference_wic)
    reference_images = reference_wic.isel(date=np.flatnonzero(reference_valid))
    reference = _robust_wic_state(reference_images, settings)

    output = {}
    for sensor, original in images.items():
        data = original.copy(deep=True)
        count = data.sizes['date']
        status = np.full(count, 'unresolved', dtype='<U16')
        reason = np.full(count, 'no_wic_support', dtype='<U32')
        status[~valid_axes[sensor]] = 'invalid_metadata'
        reason[~valid_axes[sensor]] = 'invalid_metadata'
        position_delta = np.zeros((count, 3))
        axis_delta = np.zeros((count, 3))
        axis_delta_angle = np.zeros(count)
        phase_delta = np.zeros(count)
        position_uncertainty = np.full(count, np.nan)
        axis_uncertainty = np.full(count, np.nan)
        phase_uncertainty = np.full(count, np.nan)
        coordinate_displacement = np.full(count, np.nan)
        reference_outlier = np.zeros(count, dtype=bool)
        geometry_error = np.full(count, np.nan)
        geometry_updated = np.zeros(count, dtype=bool)

        for index in np.flatnonzero(valid_axes[sensor]):
            target = (data.pointing_time.values[index].astype('datetime64[ns]').astype(np.int64) - reference['epoch'])/1e9
            state = _interpolate_wic_state(reference, target, settings['simultaneous_tolerance_seconds'])
            if state is None:
                continue
            position_uncertainty[index] = state['position_uncertainty']
            axis_uncertainty[index] = state['axis_uncertainty']
            phase_uncertainty[index] = state['phase_uncertainty']
            reference_outlier[index] = state['isolated_outlier']
            position_delta[index] = state['position'] - np.asarray(original.O_GCI.values[index], dtype=float)
            axis_delta[index] = state['axis'] - state['raw_axis']
            axis_delta_angle[index] = np.rad2deg(np.arccos(np.clip(np.dot(state['axis'], state['raw_axis']), -1, 1)))
            phase_delta[index] = state['phase'] - state['raw_phase']
            attitude_reference = np.eye(3)[np.argmin(np.abs(state['raw_axis']))]
            correction = _attitude(state['axis'], state['phase'], attitude_reference) @ _attitude(state['raw_axis'], state['raw_phase'], attitude_reference).T

            raw_dza = original.get('pointing_DZA', original.dza).values[index]
            raw_valid = np.isfinite(raw_dza) & (raw_dza >= 0) & (raw_dza < float(original.attrs.get('pointing_dzalim', 80)))
            if raw_valid.sum() < 6:
                reason[index] = 'insufficient_native_rays'
                continue
            try:
                fields, latitude, longitude, error = _correct_frame(original, index, state['position'], correction_rotation=correction)
            except ValueError:
                reason[index] = 'geometry_reconstruction_failed'
                continue
            displacement = _maximum_geographic_displacement({'glat': original.glat.values[index], 'glon': original.glon.values[index]}, fields)
            coordinate_displacement[index] = displacement
            if displacement > settings['maximum_coordinate_displacement_deg']:
                reason[index] = 'coordinate_displacement_limit'
                continue

            for name, values in fields.items():
                if name not in data:
                    data[name] = (data.sza.dims, np.full(data.sza.shape, np.nan))
                data[name].values[index] = values
            data.O_GCI.values[index] = state['position']
            data.A_GCI3.values[index] = state['axis']
            data.SPIN_PHASE.values[index] = original.SPIN_PHASE.values[index] + phase_delta[index]
            data.LAT_CTR.values[index], data.LON_CTR.values[index] = latitude, longitude
            geometry_error[index] = error
            geometry_updated[index] = True
            status[index] = 'supported'
            reason[index] = 'supported'

        for name in ['O_GCI', 'A_GCI3', 'SPIN_PHASE', 'LAT_CTR', 'LON_CTR']:
            data['pointing_original_' + name] = original[name].copy(deep=True)
        data['pointing_position_delta'] = (data.O_GCI.dims, position_delta)
        data['pointing_axis_delta'] = (data.A_GCI3.dims, axis_delta)
        data['pointing_axis_delta_angle'] = ('date', axis_delta_angle)
        data['pointing_phase_delta'] = ('date', phase_delta)
        data['pointing_position_uncertainty'] = ('date', position_uncertainty)
        data['pointing_axis_uncertainty'] = ('date', axis_uncertainty)
        data['pointing_phase_uncertainty'] = ('date', phase_uncertainty)
        data['pointing_coordinate_displacement'] = ('date', coordinate_displacement)
        data['pointing_reference_outlier'] = ('date', reference_outlier)
        data['pointing_status'] = ('date', status)
        data['pointing_status_reason'] = ('date', reason)
        data['pointing_filter_supported'] = ('date', status == 'supported')
        data['pointing_position_outlier'] = ('date', reference_outlier)
        data['pointing_geometry_updated'] = ('date', geometry_updated)
        data['pointing_baseline_angle_error'] = ('date', geometry_error)
        data.pointing_position_delta.attrs = {'units': 'km', 'long_name': 'WIC-only fitted spacecraft position minus original camera position'}
        data.pointing_axis_delta_angle.attrs = {'units': 'deg', 'long_name': 'WIC-only fitted spin-axis angular correction'}
        data.pointing_phase_delta.attrs = {'units': 'deg', 'long_name': 'WIC-only fitted spin-phase correction'}
        data.pointing_position_uncertainty.attrs = {'units': 'km'}
        data.pointing_axis_uncertainty.attrs = {'units': 'deg'}
        data.pointing_phase_uncertainty.attrs = {'units': 'deg'}
        data.pointing_coordinate_displacement.attrs = {'units': 'deg', 'long_name': 'Maximum native-pixel geographic angular displacement'}
        data.pointing_status.attrs = {'long_name': 'Robust WIC metadata pointing outcome', 'flag_meanings': 'supported unresolved invalid_metadata'}
        data.pointing_baseline_angle_error.attrs = {'units': 'deg', 'long_name': 'Maximum unchanged SZA/DZA reconstruction error'}
        for name in GEOGRAPHIC_FIELDS + MAGNETIC_FIELDS:
            data.attrs['pointing_original_' + name + '_sha256'] = hashlib.sha256(np.ascontiguousarray(original[name].values).tobytes()).hexdigest()
        data.attrs.update({
            'pointing_correction_applied': 1,
            'pointing_correction_method': 'robust_state',
            'pointing_reference_sensor': 'WIC',
            'pointing_state_model': 'constant-rate position, spin-axis and unwrapped phase; isolated leave-one-out outliers downweighted',
            'pointing_reference_matching': 'Exact-time interpolation inside supported WIC spans; endpoint tolerance only for effectively simultaneous camera frames',
            'pointing_magnetic_mapping': 'FUVVIEW3 epoch1997 lookup; original per-frame subsolar magnetic offset',
            'pointing_support': 'own pointing_valid metadata, supported WIC state, uncertainty limits and coordinate-displacement limit'
        })
        for name, value in settings.items():
            data.attrs['pointing_' + name] = value
        output[sensor] = data
    return output


def correct_pointing(images, method='clock_phase', window=11, minimum_position_outlier_km=20., reference_wic=None, state_options=None):
    """Return jointly corrected copies of a WIC/SI12/SI13 dataset mapping.

    ``method='position'`` repairs isolated spacecraft-position outliers only.
    The default ``'clock_phase'`` also smooths clock-aware spin phase with a
    local quadratic. WIC is the shared reference; camera alignments are fixed.
    The robust leave-centre-out cubic position fit and eleven-frame window
    reproduce the reviewed icBuilder metadata trial. WIC time gaps and counter
    discontinuities separate smoothing segments. Camera frames associate with
    the nearest precise WIC time, independently of their SPIN counters.
    Equidistant matches, short segments and times outside a supported WIC
    segment receive no correction or extrapolation.

    ``method='robust_state'`` is the experimental WIC-only estimator. It fits
    position/velocity, spin axis/angular rate and unwrapped phase/rate within
    observation-supported spans. Leave-one-out residuals repair isolated
    metadata errors; ambiguous transitions and insufficient support receive an
    explicit ``pointing_status='unresolved'``. ``state_options`` overrides its
    provisional numerical settings. It is opt-in while those settings are
    calibrated on the metadata census.

    ``reference_wic`` may supply the complete WIC sequence when ``images``
    contains only frames selected for the background fit. Valid spacecraft
    states from rejected WIC images still support the temporal correction;
    pixel geometry is recalculated only for the selected camera frames.

    Use ``read_idl(..., pointing_metadata=True)`` for every sensor. Counts and
    variances are unchanged. Geographic, legacy1997 Apex and background fields
    are corrected together, retaining the intersection with original support.
    Original state variables and coordinate SHA256 hashes retain provenance.
    This metadata consistency repair is not an independent attitude solution.
    """
    if not isinstance(images, Mapping) or hasattr(images, 'data_vars'):
        raise ValueError('Joint pointing correction requires a mapping of IMAGE/FUV camera datasets')
    if method not in ['position', 'clock_phase', 'robust_state']:
        raise ValueError("method must be 'position', 'clock_phase', or 'robust_state'")
    if method == 'robust_state':
        unknown = set(images) - SENSORS
        if not images or unknown:
            raise ValueError('Robust-state pointing requires one or more WIC, SI12 or SI13 camera datasets')
        if reference_wic is None and 'WIC' not in images:
            raise ValueError('Robust-state pointing requires WIC in images or as reference_wic')
        return _correct_pointing_robust_state(images, reference_wic, _state_settings(state_options))
    if set(images) != SENSORS:
        raise ValueError('Joint pointing correction requires a mapping containing WIC, SI12 and SI13')
    if state_options is not None:
        raise ValueError("state_options requires method='robust_state'")
    if not isinstance(window, int) or window < 5 or window % 2 != 1:
        raise ValueError('window must be an odd integer of at least five frames')
    if not np.isfinite(minimum_position_outlier_km) or minimum_position_outlier_km <= 0:
        raise ValueError('minimum_position_outlier_km must be finite and positive')
    valid_axes = _validate_images(images, allow_empty=reference_wic is not None)
    if reference_wic is None:
        reference_wic = images['WIC']
        reference_axes = valid_axes['WIC']
    else:
        reference_axes = _valid_wic_reference(reference_wic)
    reference_images = reference_wic.isel(date=np.flatnonzero(reference_axes))
    reference = _reference_corrections(reference_images, window, minimum_position_outlier_km)
    epoch = reference_images.pointing_time.values[0].astype('datetime64[ns]').astype(np.int64)
    output = {}
    for sensor, original in images.items():
        data = original.copy(deep=True)
        count = data.sizes['date']
        position_delta = np.zeros((count, 3))
        phase_delta = np.zeros(count)
        supported = np.zeros(count, dtype=bool)
        outlier = np.zeros(count, dtype=bool)
        geometry_error = np.full(count, np.nan)
        geometry_updated = np.zeros(count, dtype=bool)
        for index, spin in enumerate(data.SPIN.values):
            if not valid_axes[sensor][index]:
                continue
            t = (data.pointing_time.values[index].astype('datetime64[ns]').astype(np.int64) - epoch)/1e9
            right = np.searchsorted(reference['time'], t)
            candidates = np.arange(max(0, right - 1), min(len(reference['time']), right + 1))
            distances = np.abs(reference['time'][candidates] - t)
            if len(candidates) == 2 and distances[0] == distances[1]:
                continue
            ri = int(candidates[np.argmin(distances)])
            group = np.flatnonzero(reference['segment'] == reference['segment'][ri])
            if not reference['supported'][ri] or t < reference['time'][group[0]] or t > reference['time'][group[-1]]:
                continue
            supported[index] = True
            observer = np.asarray(original.O_GCI.values[index], dtype=float).copy()
            if reference['outlier'][ri]:
                coefficients = reference['polynomials'][ri]
                prediction = np.polynomial.polynomial.polyval((t - reference['time'][ri])/600, coefficients)
                observer = reference['position'][ri] + prediction
                position_delta[index] = observer - original.O_GCI.values[index]
                outlier[index] = True
            if method == 'clock_phase':
                phase_delta[index] = np.interp(t, reference['time'][group], reference['phase_delta'][group])
            if np.linalg.norm(position_delta[index]) < 1e-8 and abs(phase_delta[index]) < 1e-10:
                continue
            raw_dza = original.get('pointing_DZA', original.dza).values[index]
            raw_valid = np.isfinite(raw_dza) & (raw_dza >= 0) & (raw_dza < 75)
            if raw_valid.sum() < 6:
                # No mapped rays: a phase-only state change has no geographic
                # pixels to update. Moving the observer still needs a rotation.
                if outlier[index]:
                    raise ValueError(f'{sensor} SPIN {int(spin)} has a shared position outlier but insufficient raw coordinates to update its observer geometry')
                data.SPIN_PHASE.values[index] = original.SPIN_PHASE.values[index] + phase_delta[index]
                continue
            fields, latitude, longitude, geometry_error[index] = _correct_frame(original, index, observer, phase_delta[index])
            geometry_updated[index] = True
            for name, values in fields.items():
                if name not in data:
                    data[name] = (data.sza.dims, np.full(data.sza.shape, np.nan))
                data[name].values[index] = values
            data.O_GCI.values[index] = observer
            data.SPIN_PHASE.values[index] = original.SPIN_PHASE.values[index] + phase_delta[index]
            data.LAT_CTR.values[index], data.LON_CTR.values[index] = latitude, longitude
        for name in ['O_GCI', 'SPIN_PHASE', 'LAT_CTR', 'LON_CTR']:
            data['pointing_original_' + name] = original[name].copy(deep=True)
        data['pointing_position_delta'] = (data.O_GCI.dims, position_delta)
        data['pointing_phase_delta'] = (['date'], phase_delta)
        data['pointing_filter_supported'] = (['date'], supported)
        data['pointing_position_outlier'] = (['date'], outlier)
        data['pointing_geometry_updated'] = (['date'], geometry_updated)
        data['pointing_baseline_angle_error'] = (['date'], geometry_error)
        data.pointing_position_delta.attrs = {'units': 'km', 'long_name': 'Shared metadata spacecraft position correction'}
        data.pointing_phase_delta.attrs = {'units': 'deg', 'long_name': 'Shared clock-aware phase correction'}
        data.attrs['pointing_reference_matching'] = 'Nearest precise observation time within a supported WIC segment; equidistant matches skipped; no extrapolation'
        data.pointing_baseline_angle_error.attrs = {'units': 'deg', 'long_name': 'Maximum unchanged SZA/DZA reconstruction error'}
        for name in GEOGRAPHIC_FIELDS + MAGNETIC_FIELDS:
            data.attrs['pointing_original_' + name + '_sha256'] = hashlib.sha256(np.ascontiguousarray(original[name].values).tobytes()).hexdigest()
        data.attrs.update({'pointing_correction_applied': 1, 'pointing_correction_method': method, 'pointing_reference_sensor': 'WIC', 'pointing_window_frames': window, 'pointing_position_threshold_km': reference['threshold'], 'pointing_magnetic_mapping': 'FUVVIEW3 epoch1997 lookup; original per-frame subsolar magnetic offset', 'pointing_support': 'original mapped support intersect corrected visible shell and loader DZA limit'})
        output[sensor] = data
    return output
