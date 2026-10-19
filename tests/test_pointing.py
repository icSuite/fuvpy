from pathlib import Path
from collections import UserDict
import hashlib

import numpy as np
import pandas as pd
import pytest
from scipy.io import readsav
from scipy.spatial.transform import Rotation
import xarray as xr

from fuvpy.src.load import read_idl
from fuvpy.src.pointing import _ellipsoid, _historical_sun, _magnetic, _native_geometry, correct_pointing
from fuvpy.src.viewing_geometry import calculate_viewing_geometry
from fuvpy.utils.sunlight import subsol


def _synthetic_cameras(frames=15, spike=True, phase_noise=False):
    """Known inertial orbit and distinct native camera footprints."""
    output = {}
    for sensor, width, delay in [('WIC', 15, 0), ('SI12', 9, .01), ('SI13', 8, .02)]:
        times = pd.date_range('2001-10-21T20:00:00', periods=frames, freq='120s') + pd.to_timedelta(delay, unit='s')
        seconds = np.arange(frames)*120. + delay
        position = np.column_stack([12000 + .05*seconds + .00001*seconds**2, 2500 - .03*seconds, 33000 + .001*seconds])
        if spike:
            position[frames//2] += [300, -250, 150]
        transform = Rotation.from_euler('xyz', [.1, -.2, 40], degrees=True).as_matrix()
        glat, glon = np.meshgrid(np.linspace(60, 82, width), np.linspace(-65, 10, width), indexing='ij')
        point, _, _ = _ellipsoid(glat, glon, 130)
        inertial = point @ transform
        radial = inertial/np.linalg.norm(inertial, axis=-1, keepdims=True)
        mapped = {name: [] for name in ['sza', 'dza', 'solar_view_cosine_factor', 'solar_view_sine_factor']}
        latitude, longitude = [], []
        for time, observer in zip(times, position):
            ray = inertial - observer
            ray /= np.linalg.norm(ray, axis=-1, keepdims=True)
            mapped['sza'].append(np.rad2deg(np.arccos(np.clip(radial @ _historical_sun(time), -1, 1))))
            mapped['dza'].append(np.rad2deg(np.arccos(np.clip(np.sum(radial*(-ray), axis=-1), -1, 1))))
            geographic_observer = observer @ transform.T
            latitude.append(np.rad2deg(np.arctan2(geographic_observer[2], np.hypot(geographic_observer[0], geographic_observer[1]))))
            longitude.append(np.rad2deg(np.arctan2(geographic_observer[1], geographic_observer[0])) % 360)
            sun_latitude, sun_longitude = subsol([time.floor('s')])
            view = calculate_viewing_geometry(glat, glon, np.linalg.norm(observer), latitude[-1], longitude[-1], 130, sun_latitude[0], sun_longitude[0])
            for name in ['solar_view_cosine_factor', 'solar_view_sine_factor']:
                mapped[name].append(view[name])
        phase = .0003*seconds + (np.sin(np.arange(frames)*2)*.02 if phase_noise else 0)
        dims = ['date', 'row', 'col']
        mlat, mlon = _magnetic(glat, glon)
        fields = {name: (dims, np.broadcast_to(values, (frames, width, width)).copy()) for name, values in [('glat', glat), ('glon', glon), ('mlat', mlat), ('mlon', mlon), ('mlt', ((mlon - 20)/15) % 24)]}
        fields.update({name: (dims, np.asarray(values)) for name, values in mapped.items()})
        fields.update({'img': (dims, np.full((frames, width, width), 100., dtype=np.float32)), 'img_variance': (dims, np.full((frames, width, width), 10., dtype=np.float32)), 'O_GCI': (['date', 'pointing_xyz'], position), 'A_GCI3': (['date', 'pointing_xyz'], np.tile([0, 0, 1.], (frames, 1))), 'SPIN_PHASE': (['date'], phase), 'SPIN': (['date'], np.arange(frames) + 25000), 'LAT_CTR': (['date'], latitude), 'LON_CTR': (['date'], longitude), 'EMIS_HGT': (['date'], np.full(frames, 130.)), 'pointing_time': (['date'], times.values)})
        output[sensor] = xr.Dataset(fields, coords={'date': times.floor('s'), 'row': np.arange(width), 'col': np.arange(width), 'pointing_xyz': np.arange(3)}, attrs={'pointing_dzalim': 80})
    return output


def test_shared_position_repair_recovers_known_orbit_at_each_camera_time():
    original = _synthetic_cameras()
    frozen = {sensor: data.copy(deep=True) for sensor, data in original.items()}
    result = correct_pointing(original, method='position')
    for sensor, data in result.items():
        np.testing.assert_array_equal(data.pointing_position_outlier.values, np.arange(15) == 7)
        time = (data.pointing_time.values[7] - np.datetime64('2001-10-21T20:00:00'))/np.timedelta64(1, 's')
        expected = [12000 + .05*time + .00001*time**2, 2500 - .03*time, 33000 + .001*time]
        np.testing.assert_allclose(data.O_GCI.values[7], expected, atol=1e-9)
        assert data.pointing_baseline_angle_error.values[7] < 1e-8
        assert np.nanmax(np.abs(data.sza.values[7] - original[sensor].sza.values[7])) > .1
        np.testing.assert_array_equal(data.img, original[sensor].img)
        np.testing.assert_array_equal(data.img_variance, original[sensor].img_variance)
        xr.testing.assert_identical(original[sensor], frozen[sensor])
        # Every changed coordinate set must reproduce corrected observer geometry.
        view = calculate_viewing_geometry(data.glat.values[7], data.glon.values[7], np.linalg.norm(data.O_GCI.values[7]), data.LAT_CTR.values[7], data.LON_CTR.values[7], 130, 0, 0)
        valid = np.isfinite(data.dza.values[7])
        # The producer intersects an inertial-axis ellipsoid before its small
        # precession/nutation rotation; preserve its 0.001deg consistency check.
        np.testing.assert_allclose(view['dza'][valid], data.dza.values[7][valid], atol=1e-3)
        offset = (data.mlon.values[7][valid] - 15*data.mlt.values[7][valid]) % 360
        np.testing.assert_allclose(offset, 20, atol=1e-10)


def test_clock_phase_correction_is_shared_and_does_not_modify_inputs():
    original = _synthetic_cameras(spike=False, phase_noise=True)
    result = correct_pointing(original)
    wic = result['WIC']
    assert np.max(np.abs(wic.pointing_phase_delta.values)) > .01
    for sensor in ['SI12', 'SI13']:
        data = result[sensor]
        supported = data.pointing_filter_supported.values
        seconds = (data.pointing_time.values - wic.pointing_time.values[0])/np.timedelta64(1, 's')
        reference_seconds = (wic.pointing_time.values - wic.pointing_time.values[0])/np.timedelta64(1, 's')
        np.testing.assert_allclose(data.pointing_phase_delta.values[supported], np.interp(seconds[supported], reference_seconds, wic.pointing_phase_delta.values))
        # A sensor exposure after the last supported WIC time is not extrapolated.
        assert not supported[-1]
        np.testing.assert_array_equal(data.glat.values[-1], original[sensor].glat.values[-1])
        np.testing.assert_array_equal(data.A_GCI3, original[sensor].A_GCI3)


def test_shared_state_is_repaired_when_one_camera_has_no_science_pixels():
    source = _synthetic_cameras()
    wic = source['WIC']
    for name in ['glat', 'glon', 'sza', 'dza']:
        wic['pointing_' + name.upper()] = wic[name].copy(deep=True)
    for name in ['glat', 'glon', 'sza', 'dza', 'mlat', 'mlon', 'mlt', 'solar_view_cosine_factor', 'solar_view_sine_factor']:
        wic[name].values[7] = np.nan
    result = correct_pointing(source, method='position')
    assert result['WIC'].pointing_position_outlier.values[7]
    assert result['WIC'].pointing_geometry_updated.values[7]
    np.testing.assert_allclose(result['WIC'].O_GCI.values[7], result['SI12'].O_GCI.values[7], atol=.01)
    assert not np.any(np.isfinite(result['WIC'].sza.values[7]))


def test_short_segments_and_position_noop_keep_fields_bitwise_unchanged():
    for source in [_synthetic_cameras(frames=1), _synthetic_cameras(frames=4), _synthetic_cameras(spike=False)]:
        result = correct_pointing(source, method='position')
        for sensor, data in result.items():
            for name in source[sensor].data_vars:
                np.testing.assert_array_equal(data[name], source[sensor][name])
            assert not data.pointing_position_outlier.any()


def test_gaps_do_not_join_short_segments_or_extrapolate_position_repairs():
    source = _synthetic_cameras()
    for data in source.values():
        data.SPIN.values[8:] += 1
        data.pointing_time.values[8:] += np.timedelta64(10, 'm')
    result = correct_pointing(UserDict(source))
    for sensor, data in result.items():
        assert not data.pointing_filter_supported.any()
        np.testing.assert_array_equal(data.O_GCI, source[sensor].O_GCI)
        np.testing.assert_array_equal(data.sza, source[sensor].sza)


def test_camera_time_matching_is_independent_of_repeated_or_offset_si_counters():
    source = _synthetic_cameras(phase_noise=True)
    expected = correct_pointing(source)
    source['SI12'].SPIN.values[:] += 107
    source['SI13'].SPIN.values[8:] -= 1
    result = correct_pointing(source)
    for sensor, data in result.items():
        for name in ['pointing_filter_supported', 'pointing_position_outlier', 'pointing_phase_delta', 'O_GCI', 'glat', 'sza', 'mlat']:
            np.testing.assert_array_equal(data[name], expected[sensor][name])
        np.testing.assert_array_equal(data.SPIN, source[sensor].SPIN)


@pytest.mark.parametrize('counter_jump', [-1, -1000, 5])
def test_wic_counter_breaks_fit_each_side_separately_without_dropping_frames(counter_jump):
    source = _synthetic_cameras(frames=35, phase_noise=True)
    boundary = 13
    for data in source.values():
        data.SPIN.values[boundary:] += counter_jump
        data.SPIN_PHASE.values[boundary:] += 40
    result = correct_pointing(source)
    for section in [slice(0, boundary), slice(boundary, None)]:
        separate = correct_pointing({sensor: data.isel(date=section) for sensor, data in source.items()})
        for sensor, data in separate.items():
            joined = result[sensor].isel(date=section)
            np.testing.assert_array_equal(joined.pointing_filter_supported, data.pointing_filter_supported)
            np.testing.assert_array_equal(joined.pointing_position_outlier, data.pointing_position_outlier)
            np.testing.assert_allclose(joined.pointing_phase_delta, data.pointing_phase_delta, atol=1e-10, rtol=0)
            np.testing.assert_allclose(joined.O_GCI, data.O_GCI, atol=1e-9, rtol=0)
    for sensor, data in result.items():
        assert data.sizes['date'] == 35 and data.pointing_position_outlier.values[17]
        np.testing.assert_array_equal(data.SPIN, source[sensor].SPIN)
        np.testing.assert_array_equal(data.img, source[sensor].img)


def test_repeated_wic_counter_short_segments_keep_original_pointing():
    source = _synthetic_cameras(phase_noise=True)
    source['WIC'].SPIN.values[8:] -= 1
    result = correct_pointing(source)
    for sensor, data in result.items():
        assert not data.pointing_filter_supported.any()
        for name in source[sensor].data_vars:
            np.testing.assert_array_equal(data[name], source[sensor][name])


def test_equidistant_wic_time_match_is_left_uncorrected():
    source = _synthetic_cameras(phase_noise=True)
    time = source['WIC'].pointing_time.values
    source['SI12'].pointing_time.values[7] = time[6]+(time[7]-time[6])//2
    data = correct_pointing(source)['SI12']
    assert not data.pointing_filter_supported.values[7]
    np.testing.assert_array_equal(data.O_GCI.values[7], source['SI12'].O_GCI.values[7])
    np.testing.assert_array_equal(data.sza.values[7], source['SI12'].sza.values[7])


def test_pointing_rejects_missing_cameras_invalid_times_and_repeated_correction():
    source = _synthetic_cameras(spike=False)
    with pytest.raises(ValueError, match='WIC, SI12 and SI13'):
        correct_pointing({'WIC': source['WIC']})
    source['SI12']['id'] = 'SI13'
    with pytest.raises(ValueError, match='sensor id SI13'):
        correct_pointing(source)
    source['SI12'] = source['SI12'].drop_vars('id')
    source['SI12'].pointing_time.values[2] = np.datetime64('NaT')
    with pytest.raises(ValueError, match='finite and strictly increasing'):
        correct_pointing(source)
    result = correct_pointing(_synthetic_cameras(spike=False))
    with pytest.raises(ValueError, match='already been corrected'):
        correct_pointing(result)


def test_native_grids_reconstruct_bundled_three_sensor_coordinates():
    root = Path(__file__).parents[1]/'data/idl_files/orbit_0085'
    cameras = {}
    for folder in ['wic_data', 's12_data', 's13_data']:
        filename = next((root/folder).glob('*20001370048.idl'))
        record = readsav(filename, python_dict=True)['imageinfo'][0]
        date = pd.to_datetime(str(int(record['TIME'][0])), format='%Y%j') + pd.to_timedelta(int(record['TIME'][1]), unit='ms')
        _, transform, _, _, _, error = _native_geometry(*[np.asarray(record[name], dtype=float) for name in ['GLAT', 'GLON', 'SZA', 'DZA']], np.asarray(record['O_GCI'], dtype=float), float(record['LAT_CTR']), float(record['LON_CTR']), float(record['EMIS_HGT']), date)
        assert error < 4e-5
        np.testing.assert_allclose(transform.T @ transform, np.eye(3), atol=1e-12)
        assert np.linalg.det(transform) > .999999999
        data = read_idl(str(filename), reflat=False, pointing_metadata=True)
        cameras[str(data.id.item())] = data
    result = correct_pointing(cameras, method='position')
    for sensor, data in result.items():
        for name in cameras[sensor].data_vars:
            np.testing.assert_array_equal(data[name], cameras[sensor][name])


def test_real_shared_968_outlier_repairs_all_native_cameras():
    base = Path('/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test')
    if not (base/'wicfiles.h5').is_file():
        pytest.skip('Local private IMAGE orbit968 metadata regression corpus is unavailable')
    cameras = {}
    for sensor, prefix, middle in [('WIC', 'wic', 292), ('SI12', 's12', 291), ('SI13', 's13', 291)]:
        table = pd.read_hdf(base/f'{prefix}files.h5', key='data')
        names = table.loc[table.orbit == 968, 'filename'].iloc[middle - 7:middle + 8]
        cameras[sensor] = read_idl([str(base/f'{prefix}_data'/name) for name in names], reflat=False, viewing_geometry=True, measurement_variance=True, poisson_offset=.5, pointing_metadata=True)
    result = correct_pointing(cameras, method='position')
    for sensor, data in result.items():
        index = int(np.flatnonzero(data.SPIN.values == 25522)[0])
        assert data.pointing_position_outlier.values[index]
        assert 435 < np.linalg.norm(data.pointing_position_delta.values[index]) < 437
        if sensor == 'WIC':
            # Frozen accepted metadata prototype, not a value refit by this test.
            np.testing.assert_allclose(data.O_GCI.values[index], [11565.078542856842, 2566.268270118314, 32409.692656094223], atol=1e-8, rtol=0)
        assert data.pointing_baseline_angle_error.values[index] < 4e-5
        assert data.pointing_baseline_magnetic_error.values[index] < 1e-3
        assert np.nanmax(np.abs(data.sza.values[index] - cameras[sensor].sza.values[index])) > .1
        assert np.all(np.isfinite(data.solar_view_sine_factor.values[index][np.isfinite(data.sza.values[index])]))
        np.testing.assert_array_equal(data.img, cameras[sensor].img)
        np.testing.assert_array_equal(data.img_variance, cameras[sensor].img_variance)


def test_packaged_apex_table_retains_original_bytes():
    filename = Path(__file__).parents[1]/'fuvpy/data/pointing/mlatlon.1997a.xdr'
    assert filename.stat().st_size == 522728
    assert hashlib.sha256(filename.read_bytes()).hexdigest() == '54894228da8b60dc21cc5337cf112c4f9e40c26c65d308cf31034562f4e5b549'
