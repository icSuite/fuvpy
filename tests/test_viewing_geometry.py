from pathlib import Path
from copy import deepcopy

import numpy as np
import pytest
from scipy.io import readsav
import xarray as xr

import fuvpy.src.load as load
from fuvpy.src.load import read_idl
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


def test_unmapped_files_are_skipped_only_when_requested_and_recorded(si12_records):
    names = list(si12_records)
    with pytest.raises(ValueError, match='no valid DZA pixels'):
        read_idl(names, reflat=False, pointing_metadata=True)
    settings = {'reflat': False, 'pointing_metadata': True, 'measurement_variance': True, 'poisson_offset': .5}
    with pytest.warns(RuntimeWarning, match='Skipping unmapped SI12 frame') as messages:
        result = read_idl(names, skip_unmapped=True, **settings)
    assert len(messages) == 2
    assert 'frame_1.idl' in str(messages[0].message) and 'frame_2.idl' in str(messages[1].message)
    expected = read_idl(['frame_0.idl', 'frame_3.idl'], **settings)
    assert result.sizes['date'] == 2
    assert result.attrs.pop('skipped_unmapped_files') == 'frame_1.idl\nframe_2.idl'
    xr.testing.assert_identical(result, expected)


def test_all_unmapped_files_raise_clear_error(si12_records):
    with pytest.warns(RuntimeWarning, match='Skipping unmapped SI12 frame'):
        with pytest.raises(ValueError, match='No frames remain'):
            read_idl(['frame_1.idl', 'frame_2.idl'], reflat=False, pointing_metadata=True, skip_unmapped=True)


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
