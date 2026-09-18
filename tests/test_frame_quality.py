import numpy as np
import xarray as xr

import fuvpy.src.load as load


def _quality_dataset():
    image = np.array([
        [[np.nan]*8, [np.nan]*8],
        [[1.0, 2.0]*4, [3.0, 4.0]*4],
        [[1.0, 2.0]*4, [3.0, 4.0]*4]
    ])
    mlat = np.array([
        [[np.nan]*8, [np.nan]*8],
        [[50.0]*8, [60.0]*8],
        [[50.0]*8, [75.0]*8]
    ])
    geometry = np.where(np.isfinite(mlat), 30.0, np.nan)
    dates = np.datetime64('2001-01-01') + np.arange(3)*np.timedelta64(2, 'm')
    return xr.Dataset(
        data_vars={
            'img': (('date', 'row', 'col'), image),
            'mlat': (('date', 'row', 'col'), mlat),
            'sza': (('date', 'row', 'col'), geometry),
            'dza': (('date', 'row', 'col'), geometry),
            'hv_phos': (('date',), np.full(3, 4000.0)),
            'hv_mcp': (('date',), np.full(3, 1200.0)),
            'id': xr.DataArray('WIC')
        },
        coords={'date': dates, 'row': np.arange(2), 'col': np.arange(8)}
    )


def test_frame_quality_is_one_ordinal_variable():
    result = load._add_frame_quality(_quality_dataset())

    np.testing.assert_array_equal(result['frame_quality'], [0, 1, 2])
    np.testing.assert_array_equal(result['frame_quality'].attrs['flag_values'], [0, 1, 2])
    assert result['frame_quality'].attrs['flag_meanings'] == 'rejected usable science_ready'


def test_frame_quality_thresholds_can_be_overridden():
    result = load._add_frame_quality(_quality_dataset(), {'minimum_science_mlat': 80})

    np.testing.assert_array_equal(result['frame_quality'], [0, 1, 1])
    assert result['frame_quality'].attrs['minimum_science_mlat'] == 80


def test_detector_state_is_required_and_turn_on_pattern_is_rejected():
    base = _quality_dataset().isel(date=[2, 2, 2, 2]).copy()
    base['hv_phos'][:] = [0, 4000, 4000, 4000]
    base['hv_mcp'][:] = [0, 1200, 1200, 1200]
    base['img'][3] = np.array([[1.0, 4.0, 2.0, 3.0]*2, [4.0, 1.0, 3.0, 2.0]*2])

    result = load._add_frame_quality(base)

    np.testing.assert_array_equal(result['frame_quality'], [0, 0, 0, 2])
    assert result['turn_on_reference_correlation'][1] == 1
    assert result['turn_on_reference_correlation'][2] == 1
    assert result['turn_on_reference_correlation'][3] < 0.90
    assert result['frame_quality'].attrs['maximum_turn_on_correlation'] == 0.90

    si = base.assign(id=xr.DataArray('SI12'))
    si['hv_phos'][:] = 0
    si['hv_mcp'][:] = 4000
    result = load._add_frame_quality(si)

    np.testing.assert_array_equal(result['frame_quality'], [2, 2, 2, 2])


def test_row_spread_rejects_row_uniform_frame_without_starting_recovery():
    images = np.array([
        [[1.0]*8, [2.0]*8],
        [[0.8, 1.2]*4, [1.8, 2.2]*4],
        [[1.0, 2.0]*4, [2.0, 1.0]*4]
    ])
    mlat = np.full_like(images, 75.0)
    geometry = np.full_like(images, 30.0)
    dates = np.datetime64('2001-01-01') + np.arange(3)*np.timedelta64(2, 'm')
    data = xr.Dataset(
        data_vars={
            'img': (('date', 'row', 'col'), images),
            'mlat': (('date', 'row', 'col'), mlat),
            'sza': (('date', 'row', 'col'), geometry),
            'dza': (('date', 'row', 'col'), geometry),
            'hv_phos': (('date',), np.full(3, 4000.0)),
            'hv_mcp': (('date',), np.full(3, 1200.0)),
            'id': xr.DataArray('WIC')
        },
        coords={'date': dates, 'row': np.arange(2), 'col': np.arange(8)}
    )

    result = load._add_frame_quality(data)

    np.testing.assert_array_equal(result['frame_quality'], [0, 2, 2])
    np.testing.assert_allclose(result['median_row_relative_spread'], [0.0, 0.15, 1/3])
    assert np.all(np.isnan(result['turn_on_reference_correlation']))


def test_sparse_detector_off_frame_is_not_a_turn_on_reference():
    row_pattern = np.arange(1, 21, dtype=float)
    images = np.tile(row_pattern, (3, 20, 1))
    mlat = np.full_like(images, 75.0)
    mlat[0, 1:] = np.nan
    geometry = np.where(np.isfinite(mlat), 30.0, np.nan)
    dates = np.datetime64('2001-01-01') + np.arange(3)*np.timedelta64(2, 'm')
    data = xr.Dataset(
        data_vars={
            'img': (('date', 'row', 'col'), images),
            'mlat': (('date', 'row', 'col'), mlat),
            'sza': (('date', 'row', 'col'), geometry),
            'dza': (('date', 'row', 'col'), geometry),
            'hv_phos': (('date',), [0.0, 4000.0, 4000.0]),
            'hv_mcp': (('date',), [0.0, 1200.0, 1200.0]),
            'id': xr.DataArray('WIC')
        },
        coords={'date': dates, 'row': np.arange(20), 'col': np.arange(20)}
    )

    result = load._add_frame_quality(data)

    np.testing.assert_array_equal(result['frame_quality'], [0, 2, 2])
    assert np.isnan(result['turn_on_reference_correlation'][1])


def test_wic_reflattening_uses_quality_one_and_two(monkeypatch):
    dates = np.datetime64('2000-01-01') + np.arange(2)*np.timedelta64(2, 'm')
    images = xr.Dataset(
        data_vars={
            'img': (('date', 'row', 'col'), np.array([[[10000.0, 10000.0], [10000.0, 10000.0]], [[20.0, 20.0], [10.0, 10.0]]])),
            'sza': (('date', 'row', 'col'), np.full((2, 2, 2), 110.0)),
            'frame_quality': (('date',), np.array([0, 1], dtype=np.int8))
        },
        coords={'date': dates, 'row': np.arange(2), 'col': np.arange(2)}
    )
    monkeypatch.setattr(load, 'readsav', lambda filename: {'flatfields': np.array([[2.0, 2.0], [1.0, 1.0]])})

    result = load._reflatWIC(images)

    np.testing.assert_allclose(result['img'][1], [[10, 10], [10, 10]])
    np.testing.assert_allclose(result['img'][0], [[9990, 9990], [10000, 10000]])
