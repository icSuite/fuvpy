"""Scientific and public-interface checks for the separate two-dimensional fit."""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

import fuvpy
import fuvpy.src.background as bg
from fuvpy.src.background2d import backgroundmodel_BS2D


def _sequence(sensor='WIC', shape=(24, 20), abrupt=False):
    """Known monotonic background, independent angular curvature and bright aurora."""

    x, sine = np.meshgrid(np.linspace(-.8, 1.2, shape[0]), np.linspace(-.8, .8, shape[1]), indexing='ij')
    offsets = np.array([0, 40, -35, 35, 0]) if abrupt else np.arange(5)*8
    truth = 240 + 80*x[None] + 15*sine[None]**2 - 8*sine[None] + offsets[:, None, None]
    aurora = np.zeros_like(truth)
    aurora[:, shape[0]//3:shape[0]//3+3, shape[1]//3:shape[1]//3+2] = 500
    counts = truth + aurora + np.random.default_rng(49268).normal(0, 2, truth.shape)
    geometry = np.broadcast_to(np.rad2deg(np.arccos(x*np.cos(np.deg2rad(45)))), counts.shape)
    dates = np.datetime64('2001-01-01') + 2*np.arange(5)*np.timedelta64(1, 'm')
    images = xr.Dataset({
        'img': (('date', 'row', 'col'), counts),
        'img_variance': (('date', 'row', 'col'), np.full_like(counts, 4)),
        'sza': (('date', 'row', 'col'), geometry.copy()),
        'dza': (('date', 'row', 'col'), np.full_like(counts, 45)),
        'solar_view_sine_factor': (('date', 'row', 'col'), np.broadcast_to(sine, counts.shape).copy()),
        'id': xr.DataArray(sensor)
    }, coords={'date': dates, 'row': np.arange(shape[0]), 'col': np.arange(shape[1])})
    return images, truth, aurora


def _constant_sigma(monkeypatch, sigma=20):
    def estimate(x, counts, model, weights, knots, return_coefficients=False):
        assert return_coefficients
        return np.full_like(counts, sigma), np.array([0., sigma**2])

    monkeypatch.setattr(bg, '_noisemodel_binned', estimate)


def _fit(images, **options):
    settings = {'sKnots': [-3.5, 0, 3.5], 'dampingVal': 1e-9, 'stop': 1e-5}
    settings.update(options)
    return backgroundmodel_BS2D(images, **settings)


@pytest.mark.parametrize('freeze_iteration', [1, 2, 3])
def test_requested_sigma_estimate_is_frozen_and_consumed(monkeypatch, freeze_iteration):
    estimates = []

    def estimate(x, counts, model, weights, knots, return_coefficients=False):
        assert return_coefficients
        sigma = 30 - 5*len(estimates)
        estimates.append({'weights': weights.copy(), 'sigma': sigma})
        return np.full_like(counts, sigma), np.array([0., sigma**2])

    monkeypatch.setattr(bg, '_noisemodel_binned', estimate)
    images, unused, aurora = _sequence()
    result = _fit(images, dgsigma_freeze_iteration=freeze_iteration, stop=1e8)

    assert len(estimates) == freeze_iteration
    assert result.dgfit_iterations.item() == freeze_iteration + 1
    if freeze_iteration > 1:
        assert np.any(estimates[1]['weights'] < estimates[0]['weights'])
    np.testing.assert_array_equal(result.dgsigma, estimates[-1]['sigma'])
    np.testing.assert_array_equal(result.dgsigma_independent, estimates[-1]['sigma'])
    np.testing.assert_array_equal(result.dgmodel.attrs['sigma_model_coefficients'], [0, estimates[-1]['sigma']**2])
    assert result.dgmodel.attrs['dgsigma_freeze_iteration'] == freeze_iteration
    assert result.dgsigma_fit_coefficients.shape == (freeze_iteration, 2)
    for suffix in ['', '_independent']:
        expected = (1-np.minimum(1, (result['dgimg'+suffix]/(5*result['dgsigma'+suffix]))**2))**2
        np.testing.assert_allclose(result['dgweight'+suffix], expected)


def test_independent_recovers_abrupt_background_and_preserves_aurora(monkeypatch):
    _constant_sigma(monkeypatch)
    images, truth, aurora = _sequence(abrupt=True)
    original = images.copy(deep=True)
    result = _fit(images)

    independent_error = np.sqrt(np.mean((result.dgmodel_independent.values-truth)**2))
    temporal_error = np.sqrt(np.mean((result.dgmodel.values-truth)**2))
    assert independent_error < 1
    assert temporal_error > 10*independent_error
    assert np.median(result.dgimg_independent.values[aurora > 0]) > 490
    assert result.dgfit_converged.item()
    assert result.dgfit_converged_independent.values.all()
    np.testing.assert_allclose(result.dgimg, images.img-result.dgmodel)
    np.testing.assert_allclose(result.dgimg_independent, images.img-result.dgmodel_independent)
    assert np.max(np.abs(result.dgmodel.values-result.dgmodel_independent.values)) > 20
    assert np.all(np.diff(result.dgmodel_coefficients.values, axis=0) >= -1e-7)
    assert np.all(np.diff(result.dgmodel_independent_coefficients.values, axis=1) >= -1e-7)
    xr.testing.assert_identical(images, original)


def test_temporal_surface_recovers_known_angular_and_time_variation(monkeypatch):
    _constant_sigma(monkeypatch)
    images, truth, aurora = _sequence()
    result = _fit(images)
    assert np.sqrt(np.mean((result.dgmodel.values-truth)**2)) < 1
    assert np.median(result.dgimg.values[aurora > 0]) > 490


def test_inherited_background_and_weights_are_ignored(monkeypatch):
    _constant_sigma(monkeypatch)
    images, unused, aurora = _sequence()
    reference = _fit(images)
    contaminated = images.copy(deep=True)
    for suffix in ['', '_independent']:
        contaminated['dgmodel'+suffix] = xr.full_like(images.img, -1e8)
        contaminated['dgweight'+suffix] = xr.zeros_like(images.img)
        contaminated['dgsigma'+suffix] = xr.full_like(images.img, 1e8)
    repeated = _fit(contaminated)
    for suffix in ['', '_independent']:
        for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma']:
            np.testing.assert_array_equal(repeated[name+suffix], reference[name+suffix])


def test_rerun_with_new_knots_and_sigma_freeze_replaces_old_fit_state(monkeypatch):
    _constant_sigma(monkeypatch)
    images, unused, aurora = _sequence()
    first = _fit(images)
    options = {'sKnots': [-3.5, -.5, 0, .5, 3.5], 'dgsigma_freeze_iteration': 3, 'temporal_knot_spacing': 3}
    rerun = _fit(first, **options)
    clean = _fit(images, **options)
    for name in ['dgmodel', 'dgmodel_independent', 'dgmodel_coefficients', 'dgmodel_independent_coefficients', 'dgsigma_fit_coefficients']:
        xr.testing.assert_identical(rerun[name], clean[name])
    assert first.dgsigma_fit_coefficients.shape == (2, 2)
    assert rerun.dgsigma_fit_coefficients.shape == (3, 2)


def test_independent_sigma_uses_own_prediction_not_temporal_map(monkeypatch):
    def estimate(x, counts, model, weights, knots, return_coefficients=False):
        return np.sqrt(.003*model**2+16), np.array([.003, 16.])

    monkeypatch.setattr(bg, '_noisemodel_binned', estimate)
    images, unused, aurora = _sequence(abrupt=True)
    result = _fit(images)
    np.testing.assert_allclose(result.dgsigma_independent**2, .003*result.dgmodel_independent**2+16)
    assert np.max(np.abs(result.dgsigma_independent.values-result.dgsigma.values)) > .1


@pytest.mark.parametrize('sensor', ['SI12', 'SI13'])
def test_native_si_shape_quality_and_unsupported_pixels(monkeypatch, sensor):
    _constant_sigma(monkeypatch)
    images, unused, aurora = _sequence(sensor, shape=(128, 128))
    images['frame_quality'] = ('date', [1, 2, 2, 2, 2])
    images.img.values[:, 0, 0] = np.nan
    images.sza.values[:, 0, 1] = -1
    images.solar_view_sine_factor.values[:, 0, 2] = 1.1
    result = _fit(images)
    assert result.img.shape == (4, 128, 128)
    np.testing.assert_array_equal(result.date, images.date.values[1:])
    for suffix in ['', '_independent']:
        for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma']:
            assert np.isnan(result[name+suffix].values[:, 0, :3]).all()
            assert np.isfinite(result[name+suffix].values[:, 1:, :]).all()


def test_single_image_retains_its_observation_time(monkeypatch):
    _constant_sigma(monkeypatch)
    images, unused, aurora = _sequence()
    frame = images.isel(date=2)
    result = _fit(frame)
    assert result.sizes['date'] == 1
    assert result.date.values[0] == images.date.values[2]
    np.testing.assert_allclose(result.dgmodel, result.dgmodel_independent, atol=1e-8)


@pytest.mark.parametrize('options, message', [
    ({'dgsigma_freeze_iteration': 0}, 'dgsigma_freeze_iteration'),
    ({'dgsigma_freeze_iteration': True}, 'dgsigma_freeze_iteration'),
    ({'dgsigma_freeze_iteration': 1.5}, 'dgsigma_freeze_iteration'),
    ({'sKnots': [0, 0, 1]}, 'sKnots'),
    ({'sKnots': [0, np.nan, 1]}, 'sKnots'),
    ({'sine_knots': [-1, 0, 0, 1]}, 'sine_knots'),
    ({'sine_knots': [1, -1]}, 'sine_knots'),
    ({'temporal_knot_spacing': 0}, 'temporal.knot.spacing'),
    ({'variance_weighted': True}, 'img_variance')
])
def test_invalid_controls_are_rejected(options, message):
    images, unused, aurora = _sequence()
    if message == 'img_variance':
        images = images.drop_vars('img_variance')
    with pytest.raises((ValueError, TypeError), match=message):
        _fit(images, **options)


def test_accessor_export_and_both_results_survive_netcdf(monkeypatch, tmp_path):
    _constant_sigma(monkeypatch)
    images, unused, aurora = _sequence()
    assert fuvpy.backgroundmodel_BS2D is backgroundmodel_BS2D
    fitted = images.fuv.backgroundmodel_BS2D(sKnots=[-3.5, 0, 3.5], dampingVal=1e-9)
    path = tmp_path/'both_background_models.nc'
    fitted.to_netcdf(path)
    with xr.open_dataset(path) as restored:
        for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma', 'dgmodel_independent', 'dgimg_independent', 'dgweight_independent', 'dgsigma_independent', 'dgmodel_coefficients', 'dgmodel_independent_coefficients', 'dgsigma_fit_coefficients', 'dgfit_iterations_independent', 'dgfit_converged_independent']:
            np.testing.assert_array_equal(restored[name], fitted[name])
        np.testing.assert_array_equal(restored.dgmodel.attrs['sigma_model_coefficients'], [0, 400])
        assert restored.dgmodel.attrs['dgsigma_freeze_iteration'] == 2


def test_original_1d_defaults_match_frozen_preimplementation_result():
    """Fixture uses git65603e0 background.py, _synthetic_wic and no fit options.

    The archived implementation was executed in an independent namespace;
    the fixture therefore does not compare two calls to the modified code.
    """
    from test_background_dgsigma_freeze import _synthetic_wic

    fitted = bg.backgroundmodel_BS(_synthetic_wic())
    archive = Path(__file__).parent/'data'/'background_1d_default_65603e0.npz'
    with np.load(archive) as baseline:
        assert baseline['source_commit'].item() == '65603e0df028ea41e9cb607a0310f682fa641954'
        for name in ['dgmodel', 'dgimg', 'dgsigma', 'dgweight', 'dgnorm_residual', 'dgnorm_model', 'dgsolver_flag', 'dgfit_iterations', 'dgfit_converged']:
            np.testing.assert_array_equal(fitted[name].values, baseline[name])


def test_both_public_methods_apply_one_joint_pointing_correction(monkeypatch, tmp_path):
    import fuvpy.src.pointing as pointing
    from test_pointing import _synthetic_cameras

    source = _synthetic_cameras()
    for sensor, data in source.items():
        data['id'] = sensor
        data['source_counts'] = data.img.copy(deep=True)*2
    original = {sensor: data.copy(deep=True) for sensor, data in source.items()}
    calls = []
    real_corrector = pointing.correct_pointing

    def record(images, **options):
        result = real_corrector(images, **options)
        calls.append(result)
        return result

    def constant_sigma(x, counts, model, weights, knots, return_coefficients=False):
        sigma = np.full_like(counts, 20)
        return (sigma, np.array([0., 400.])) if return_coefficients else sigma

    monkeypatch.setattr(pointing, 'correct_pointing', record)
    monkeypatch.setattr(bg, '_noisemodel_binned', constant_sigma)
    results = []
    for method in [fuvpy.backgroundmodel_BS, fuvpy.backgroundmodel_BS2D]:
        fitted = method(source, pointing_correction=True, pointing_options={'method': 'position'}, sKnots=[-3.5, 0, 3.5], variance_weighted=True)
        results.append(fitted)
        assert len(calls) == len(results)
        for sensor, data in fitted.items():
            assert data.pointing_position_outlier.values[7]
            assert np.linalg.norm(data.O_GCI.values[7]-source[sensor].O_GCI.values[7]) > 400
            for name in ['img', 'source_counts', 'img_variance']:
                np.testing.assert_array_equal(data[name], source[sensor][name])
            assert 'dgmodel' in data
            assert ('dgmodel_independent' in data) == (method is fuvpy.backgroundmodel_BS2D)
            for name in ['sza', 'dza', 'glat', 'glon', 'solar_view_sine_factor']:
                np.testing.assert_array_equal(data[name], calls[-1][sensor][name])
            xr.testing.assert_identical(source[sensor], original[sensor])
    for sensor in source:
        for name in ['O_GCI', 'sza', 'dza', 'glat', 'glon', 'solar_view_sine_factor']:
            np.testing.assert_array_equal(results[0][sensor][name], results[1][sensor][name])
    filename = tmp_path/'joint_pointing_and_both_backgrounds.nc'
    results[1]['WIC'].to_netcdf(filename)
    with xr.open_dataset(filename) as restored:
        assert restored.attrs['pointing_correction_applied'] == 1
        for name in ['O_GCI', 'pointing_position_delta', 'dgmodel', 'dgmodel_independent', 'img_variance', 'source_counts']:
            np.testing.assert_array_equal(restored[name], results[1]['WIC'][name])


@pytest.mark.parametrize('method', [fuvpy.backgroundmodel_BS, fuvpy.backgroundmodel_BS2D])
def test_pointing_requires_joint_mapping_but_ordinary_mapping_needs_no_metadata(monkeypatch, method):
    def constant_sigma(x, counts, model, weights, knots, return_coefficients=False):
        sigma = np.full_like(counts, 20)
        return (sigma, np.array([0., 400.])) if return_coefficients else sigma

    monkeypatch.setattr(bg, '_noisemodel_binned', constant_sigma)
    images, unused, aurora = _sequence()
    with pytest.raises(ValueError, match='mapping'):
        method(images, pointing_correction=True)
    sensors = {sensor: _sequence(sensor)[0] for sensor in ['WIC', 'SI12', 'SI13']}
    fitted = method(sensors, sKnots=[-3.5, 0, 3.5])
    assert set(fitted) == set(sensors)
    assert all('dgmodel' in data for data in fitted.values())
    assert all(not data.attrs.get('pointing_correction_applied') for data in fitted.values())


def test_nonconvergence_raises_without_modifying_input(monkeypatch):
    _constant_sigma(monkeypatch)
    images, unused, aurora = _sequence()
    original = images.copy(deep=True)
    with pytest.raises(RuntimeError, match='did not converge'):
        _fit(images, max_iterations=2, dgsigma_freeze_iteration=1, stop=1e-15)
    xr.testing.assert_identical(images, original)


def test_all_rejected_pixels_raise_without_returning_ridge_only_solution(monkeypatch):
    _constant_sigma(monkeypatch, sigma=1e-12)
    images, unused, aurora = _sequence()
    original = images.copy(deep=True)
    with pytest.raises(RuntimeError, match='rejected all supported pixels'):
        _fit(images, dampingVal=100)
    xr.testing.assert_identical(images, original)
