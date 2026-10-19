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


def _validate_images(images):
    required = ['pointing_time', 'O_GCI', 'A_GCI3', 'SPIN_PHASE', 'SPIN', 'LAT_CTR', 'LON_CTR', 'EMIS_HGT'] + GEOGRAPHIC_FIELDS + MAGNETIC_FIELDS + SOLAR_FIELDS
    for sensor, data in images.items():
        if 'id' in data and str(data.id.item()) != sensor:
            raise ValueError(f'{sensor} mapping entry contains sensor id {data.id.item()}')
        missing = [name for name in required if name not in data]
        if missing:
            raise ValueError(f'{sensor} pointing metadata missing: {", ".join(missing)}; load with pointing_metadata=True')
        if data.sizes.get('date', 0) < 1:
            raise ValueError(f'{sensor} requires at least one frame')
        time = data.pointing_time.values.astype('datetime64[ns]')
        if np.any(np.isnat(time)) or np.any(np.diff(time.astype(np.int64)) <= 0):
            raise ValueError(f'{sensor} pointing_time must be finite and strictly increasing')
        for name in required[1:8]:
            if not np.all(np.isfinite(data[name].values)):
                raise ValueError(f'{sensor} {name} metadata must be finite')
        spin = data.SPIN.values
        if np.any(spin != np.rint(spin)):
            raise ValueError(f'{sensor} SPIN must contain integers')
        norm = np.linalg.norm(data.A_GCI3.values, axis=1)
        if np.any(np.abs(norm - 1) > 1e-3):
            raise ValueError(f'{sensor} A_GCI3 must contain unit spin-axis vectors')
        if data.attrs.get('pointing_correction_applied', False):
            raise ValueError(f'{sensor} pointing has already been corrected')


def _correct_frame(data, index, observer, phase_delta):
    raw = {name: np.asarray(data.get('pointing_' + name.upper(), data[name]).isel(date=index).values, dtype=float) for name in GEOGRAPHIC_FIELDS}
    height = float(data.EMIS_HGT.values[index])
    time = data.pointing_time.values[index]
    original_observer = np.asarray(data.O_GCI.values[index], dtype=float)
    look, transform, sun, equator, ratio, error = _native_geometry(raw['glat'], raw['glon'], raw['sza'], raw['dza'], original_observer, float(data.LAT_CTR.values[index]), float(data.LON_CTR.values[index]), height, time)
    axis = np.asarray(data.A_GCI3.values[index], dtype=float).copy()
    axis /= np.linalg.norm(axis)
    look = look @ Rotation.from_rotvec(np.deg2rad(phase_delta)*axis).as_matrix().T
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
    # Apex longitude is ill-conditioned at the magnetic pole; compare unit vectors.
    old_mlat, old_mlon = _magnetic(raw['glat'], raw['glon'])
    baseline = np.isfinite(data.mlat.values[index]) & np.isfinite(data.mlon.values[index]) & np.isfinite(data.mlt.values[index])
    def unit(latitude, longitude):
        latitude, longitude = np.deg2rad([latitude, longitude])
        return np.stack([np.cos(latitude)*np.cos(longitude), np.cos(latitude)*np.sin(longitude), np.sin(latitude)], axis=-1)
    magnetic_error = np.rad2deg(np.linalg.norm(unit(old_mlat[baseline], old_mlon[baseline]) - unit(data.mlat.values[index][baseline], data.mlon.values[index][baseline]), axis=-1)).max() if np.any(baseline) else np.nan
    if magnetic_error > 1e-3:
        raise ValueError(f'Stored magnetic coordinates disagree with the historical Apex grid by {magnetic_error:.6g} degrees')
    offset = (data.mlon.values[index][baseline] - 15*data.mlt.values[index][baseline] + 180) % 360 - 180
    centre_offset = np.rad2deg(np.angle(np.mean(np.exp(1j*np.deg2rad(offset))))) if len(offset) else np.nan
    if len(offset) and np.max(np.abs((offset - centre_offset + 180) % 360 - 180)) > 1e-3:
        raise ValueError('Stored magnetic longitude and local time do not share one subsolar offset')
    mlat, mlon = _magnetic(glat, glon)
    fields = {'glat': glat, 'glon': glon, 'sza': sza, 'dza': dza, 'mlat': mlat, 'mlon': mlon, 'mlt': ((mlon - centre_offset)/15) % 24}
    fields.update({name: view[name] for name in SOLAR_FIELDS})
    fields = {name: np.where(support, values, np.nan) for name, values in fields.items()}
    return fields, centre_latitude, centre_longitude, error, float(magnetic_error)


def correct_pointing(images, method='clock_phase', window=11, minimum_position_outlier_km=20.):
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

    Use ``read_idl(..., pointing_metadata=True)`` for every sensor. Counts and
    variances are unchanged. Geographic, legacy1997 Apex and background fields
    are corrected together, retaining the intersection with original support.
    Original state variables and coordinate SHA256 hashes retain provenance.
    This metadata consistency repair is not an independent attitude solution.
    """
    if not isinstance(images, Mapping) or hasattr(images, 'data_vars') or set(images) != SENSORS:
        raise ValueError('Joint pointing correction requires a mapping containing WIC, SI12 and SI13')
    if method not in ['position', 'clock_phase']:
        raise ValueError("method must be 'position' or 'clock_phase'")
    if not isinstance(window, int) or window < 5 or window % 2 != 1:
        raise ValueError('window must be an odd integer of at least five frames')
    if not np.isfinite(minimum_position_outlier_km) or minimum_position_outlier_km <= 0:
        raise ValueError('minimum_position_outlier_km must be finite and positive')
    _validate_images(images)
    reference = _reference_corrections(images['WIC'], window, minimum_position_outlier_km)
    epoch = images['WIC'].pointing_time.values[0].astype('datetime64[ns]').astype(np.int64)
    output = {}
    for sensor, original in images.items():
        data = original.copy(deep=True)
        count = data.sizes['date']
        position_delta = np.zeros((count, 3))
        phase_delta = np.zeros(count)
        supported = np.zeros(count, dtype=bool)
        outlier = np.zeros(count, dtype=bool)
        geometry_error = np.full(count, np.nan)
        magnetic_error = np.full(count, np.nan)
        geometry_updated = np.zeros(count, dtype=bool)
        for index, spin in enumerate(data.SPIN.values):
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
            fields, latitude, longitude, geometry_error[index], magnetic_error[index] = _correct_frame(original, index, observer, phase_delta[index])
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
        data['pointing_baseline_magnetic_error'] = (['date'], magnetic_error)
        data.pointing_position_delta.attrs = {'units': 'km', 'long_name': 'Shared metadata spacecraft position correction'}
        data.pointing_phase_delta.attrs = {'units': 'deg', 'long_name': 'Shared clock-aware phase correction'}
        data.attrs['pointing_reference_matching'] = 'Nearest precise observation time within a supported WIC segment; equidistant matches skipped; no extrapolation'
        data.pointing_baseline_angle_error.attrs = {'units': 'deg', 'long_name': 'Maximum unchanged SZA/DZA reconstruction error'}
        data.pointing_baseline_magnetic_error.attrs = {'units': 'deg', 'long_name': 'Maximum unchanged magnetic unit-vector direction error'}
        for name in GEOGRAPHIC_FIELDS + MAGNETIC_FIELDS:
            data.attrs['pointing_original_' + name + '_sha256'] = hashlib.sha256(np.ascontiguousarray(original[name].values).tobytes()).hexdigest()
        data.attrs.update({'pointing_correction_applied': 1, 'pointing_correction_method': method, 'pointing_reference_sensor': 'WIC', 'pointing_window_frames': window, 'pointing_position_threshold_km': reference['threshold'], 'pointing_magnetic_mapping': 'FUVVIEW3 epoch1997 lookup; original per-frame subsolar magnetic offset', 'pointing_support': 'original mapped support intersect corrected visible shell and loader DZA limit'})
        output[sensor] = data
    return output
