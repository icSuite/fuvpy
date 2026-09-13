import numpy as np
import pytest
import xarray as xr

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic


def _directional_background(spencer_taper=False):
    """Return a small background that follows the retained model exactly."""

    n_time = 3
    n_row = 15
    n_col = 15

    sza_axis = np.linspace(20, 150, n_row)
    dza_axis = np.linspace(5, 75, n_col)
    sza_0, dza_0 = np.meshgrid(sza_axis, dza_axis, indexing='ij')
    sza = np.stack([sza_0]*n_time)
    dza = np.stack([dza_0]*n_time)
    phase = np.linspace(-np.pi, np.pi, n_col)[None, None, :]
    sine_factor = np.sin(np.deg2rad(sza))*np.sin(np.deg2rad(dza))*np.sin(phase)

    fraction = np.cos(np.deg2rad(sza))/np.cos(np.deg2rad(dza))
    taper = background._spencer_taper_weights(sza) if spencer_taper else 1
    image = 100 + 20*fraction + taper*sine_factor*(5 + 2*fraction)
    dates = np.datetime64('2001-01-01') + 2*np.arange(n_time)*np.timedelta64(1, 'm')

    return xr.Dataset(
        data_vars={
            'img': (('date', 'row', 'col'), image),
            'sza': (('date', 'row', 'col'), sza),
            'dza': (('date', 'row', 'col'), dza),
            'solar_view_sine_factor': (('date', 'row', 'col'), sine_factor),
            'solar_view_cosine_factor': (('date', 'row', 'col'), np.zeros_like(sine_factor)),
            'id': xr.DataArray('WIC')
        },
        coords={'date': dates, 'row': np.arange(n_row), 'col': np.arange(n_col)}
    )


def _fit_exact_background(monkeypatch, images, **options):
    monkeypatch.setattr(background, '_noisemodel_binned', lambda fraction, image, model, weights, knots: np.full_like(image, 1e6))
    return background.backgroundmodel_BS(images, spatial_model='legacy_x_azimuth', sKnots=[-3.5, 0, 3.5], sOrder=1, tukeyVal=1e6, stop=1e-8, **options)


def test_directional_model_recovers_background_with_and_without_taper(monkeypatch):
    for spencer_taper in [False, True]:
        images = _directional_background(spencer_taper)
        result = _fit_exact_background(monkeypatch, images, spencer_taper=spencer_taper)

        assert np.sqrt(np.nanmean(result['dgimg'].values**2)) < 1e-8
        np.testing.assert_allclose(result['dgmodel'], result['dgmodel_base'] + result['dgmodel_sine'])
        assert result['dgmodel'].attrs['spencer_taper'] == ('enabled' if spencer_taper else 'disabled')
        assert ('spencer_taper_factor' in result) == spencer_taper


def test_directional_knots_default_to_base_knots(monkeypatch):
    result = _fit_exact_background(monkeypatch, _directional_background())

    np.testing.assert_array_equal(result['dgmodel'].attrs['directional_x_knots'], [-3.5, 0, 3.5])


def test_directional_knots_must_be_increasing_and_finite():
    for knots in [[0, 0, 1], [1, 0], [0, np.nan, 1], [[0, 1]]]:
        with pytest.raises(ValueError, match='directional_x_knots'):
            background.backgroundmodel_BS(_directional_background(), spatial_model='legacy_x_azimuth', directional_x_knots=knots)


def test_directional_options_and_geometry_are_validated():
    images = _synthetic_wic()

    with pytest.raises(ValueError, match='directional_x_knots'):
        background.backgroundmodel_BS(images, directional_x_knots=[-3.5, 0, 3.5])
    with pytest.raises(ValueError, match='spencer_taper'):
        background.backgroundmodel_BS(images, spencer_taper=True)
    with pytest.raises(ValueError, match='viewing_geometry=True'):
        background.backgroundmodel_BS(images, spatial_model='legacy_x_azimuth')


def test_monotonic_constraint_applies_to_base(monkeypatch):
    images = _directional_background()
    result = _fit_exact_background(monkeypatch, images, monotonic=True)
    fraction = np.cos(np.deg2rad(images['sza']))/np.cos(np.deg2rad(images['dza']))

    for frame in range(images.sizes['date']):
        order = np.argsort(fraction.values[frame].ravel())
        ordered_base = result['dgmodel_base'].values[frame].ravel()[order]
        ordered_base = ordered_base[np.isfinite(ordered_base)]
        assert np.all(np.diff(ordered_base) >= -1e-8)

    assert result['dgmodel'].attrs['monotonic_constraint'] == 'legacy_base_nondecreasing'


def test_directional_outputs_survive_netcdf_round_trip(monkeypatch, tmp_path):
    images = _directional_background(spencer_taper=True)
    result = _fit_exact_background(monkeypatch, images, spencer_taper=True)
    filename = tmp_path/'directional.nc'

    result.to_netcdf(filename)
    loaded = xr.load_dataset(filename)

    for field in ['solar_view_cosine_factor', 'solar_view_sine_factor', 'spencer_taper_factor', 'dgmodel_base', 'dgmodel_sine']:
        np.testing.assert_allclose(loaded[field], result[field])
    assert loaded['dgmodel'].attrs['spatial_model'] == 'legacy_x_azimuth'
    assert loaded['dgmodel'].attrs['spencer_taper'] == 'enabled'
    np.testing.assert_array_equal(loaded['dgmodel'].attrs['spatial_knots'], [-3.5, 0, 3.5])
    np.testing.assert_array_equal(loaded['dgmodel'].attrs['directional_x_knots'], [-3.5, 0, 3.5])
    np.testing.assert_array_equal(loaded['dgmodel'].attrs['temporal_knots_minutes_from_first_image'], [0, 4])
    assert loaded['dgmodel'].attrs['spatial_order'] == 1
    assert loaded['dgmodel'].attrs['fit_iterations'] == result['dgmodel'].attrs['fit_iterations']
    for name in ['temporal_knot_count', 'temporal_order', 'noise_model', 'spatial_weighting', 'dgsigma_estimation', 'robust_weighting', 'tukey_value', 'measurement_variance', 'variance_weighted', 'monotonic_requested', 'monotonic_constraint', 'monotonic_solver', 'monotonic_fallback', 'damping_value', 'damping_scale', 'convergence_tolerance', 'maximum_iterations', 'fit_converged']:
        assert loaded['dgmodel'].attrs[name] == result['dgmodel'].attrs[name]
    assert int(loaded['dgsolver_flag']) == int(result['dgsolver_flag'])
    np.testing.assert_array_equal(loaded['dgsolver_flag'].attrs['flag_values'], [0, 1, 2, 3, 4])
