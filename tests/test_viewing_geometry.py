from pathlib import Path
from copy import deepcopy

import numpy as np
import pytest
from scipy.io import readsav
import xarray as xr

import fuvpy.src.load as load
from fuvpy.src.load import _pointing_record_valid, read_idl
from fuvpy.src.viewing_geometry import calculate_viewing_geometry


def test_solar_viewing_factors_are_independent_of_polar_longitude():
    geometry = calculate_viewing_geometry([90, 90], [0, 180], 50000, 60, 20, 130, 20, 30)

    np.testing.assert_allclose(geometry['solar_view_cosine_factor'][0], geometry['solar_view_cosine_factor'][1], atol=1e-14)
    np.testing.assert_allclose(geometry['solar_view_sine_factor'][0], geometry['solar_view_sine_factor'][1], atol=1e-14)


def test_bundled_idl_geometry_reproduces_stored_dza():
    for name in ['wic_data', 's12_data', 's13_data']:
        directory = Path(__file__).parents[1]/'data/idl_files/orbit_0085'/name
        filename = next(directory.glob('*20001370048.idl'))
        image = read_idl(str(filename), reflat=False, viewing_geometry=True)
        valid = np.isfinite(image['dza'].values)

        assert np.any(valid)
        assert np.all(np.isfinite(image['solar_view_cosine_factor'].values[valid]))
        assert np.all(np.isfinite(image['solar_view_sine_factor'].values[valid]))
        assert image['solar_view_cosine_factor'].shape == image['img'].shape
        assert image['solar_view_sine_factor'].shape == image['img'].shape


@pytest.fixture
def si12_records(monkeypatch):
    directory = Path(__file__).parents[1]/'data/idl_files/orbit_0085/s12_data'
    original = readsav(str(next(directory.glob('*20001370048.idl'))))['imageinfo']
    records = {f'frame_{index}.idl': deepcopy(original) for index in range(4)}
    for index, record in enumerate(records.values()):
        record['time'][0][1] += index*120000
    for name in ['frame_1.idl', 'frame_2.idl']:
        records[name]['dza'][0][:] = -1
        records[name]['glat'][0][:] = -1e31
        records[name]['glon'][0][:] = -1e31
    monkeypatch.setattr(load, 'readsav', lambda name: {'imageinfo': records[name]})
    return records


def test_unmapped_files_remain_available_as_metadata_only_frames(si12_records):
    names = list(si12_records)
    with pytest.raises(ValueError, match='no valid DZA pixels'):
        read_idl(names, reflat=False, pointing_metadata=True)
    settings = {'reflat': False, 'pointing_metadata': True, 'measurement_variance': True, 'poisson_offset': .5}
    result = read_idl(names, skip_unmapped=True, **settings)
    assert result.sizes['date'] == 4
    assert result.attrs['metadata_only_files'] == 'frame_1.idl\nframe_2.idl'
    np.testing.assert_array_equal(result.pointing_valid, True)
    np.testing.assert_array_equal(result.frame_quality.values[1:3], 0)
    assert not np.any(np.isfinite(result.dza.values[1:3]))


def test_all_unmapped_files_remain_metadata_only(si12_records):
    result = read_idl(['frame_1.idl', 'frame_2.idl'], reflat=False, pointing_metadata=True, skip_unmapped=True)
    np.testing.assert_array_equal(result.pointing_valid, True)
    np.testing.assert_array_equal(result.frame_quality, 0)


@pytest.mark.parametrize('problem', ['axis', 'phase', 'position', 'latitude', 'height', 'alignment'])
def test_corrupt_pointing_fields_are_flagged_per_frame(si12_records, problem):
    record = si12_records['frame_0.idl']
    if problem == 'axis':
        record['a_gci3'][0][:] = .7247730218805373
    elif problem == 'phase':
        record['spin_phase'][0] = np.nan
    elif problem == 'position':
        record['o_gci'][0][:] = 0
    elif problem == 'latitude':
        record['lat_ctr'][0] = 95
    elif problem == 'height':
        record['emis_hgt'][0] = np.nan
    else:
        record['inst_roll'][0] = np.nan
    result = read_idl(['frame_0.idl', 'frame_3.idl'], reflat=False, pointing_metadata=True, skip_unmapped=True)
    np.testing.assert_array_equal(result.pointing_valid.values, [False, True])
    assert result.frame_quality.values[0] == 2


def test_fractional_spin_is_invalid_but_repeated_integer_spin_is_valid(si12_records):
    result = read_idl(['frame_0.idl', 'frame_3.idl'], reflat=False, pointing_metadata=True)
    assert result.SPIN.values[0] == result.SPIN.values[1]
    np.testing.assert_array_equal(result.pointing_valid, True)
    result.SPIN.values[0] += .5
    assert not _pointing_record_valid(result.isel(date=[0]))


def test_unreadable_pointing_time_is_skipped_with_filename(si12_records):
    si12_records['frame_0.idl']['time'][0][1] = 100000000
    with pytest.warns(RuntimeWarning, match='invalid TIME'):
        result = read_idl(['frame_0.idl', 'frame_3.idl'], reflat=False, pointing_metadata=True, skip_unmapped=True)
    assert result.sizes['date'] == 1
    assert result.attrs['skipped_invalid_time_files'] == 'frame_0.idl'
    assert result.attrs['pointing_invalid_time_audit'].startswith('frame_0.idl\tpointing_valid=false\treason=invalid_time\t')


@pytest.mark.parametrize('problem', ['mismatch', 'stored_only', 'reconstructed_only'])
def test_skip_unmapped_does_not_hide_partial_or_inconsistent_geometry(si12_records, problem):
    record = si12_records['frame_0.idl']
    if problem == 'mismatch':
        record['dza'][0][record['dza'][0] >= 0] += 1
    elif problem == 'stored_only':
        record['glat'][0][:] = -1e31
        record['glon'][0][:] = -1e31
    else:
        record['dza'][0][:] = -1
    with pytest.raises(ValueError, match='Viewing geometry could not|Reconstructed viewing geometry disagrees'):
        read_idl(['frame_0.idl'], reflat=False, pointing_metadata=True, skip_unmapped=True)
