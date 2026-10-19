import numpy as np
import pytest
import xarray as xr

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic
from test_background_directional import _directional_background


def closeup(aurora=0, offset=0):
    x = np.broadcast_to(np.linspace(.08, .22, 40), (20, 40))
    truth = 200 + 80*x + offset
    sza = np.rad2deg(np.arccos(x*np.cos(np.deg2rad(45))))
    image = xr.Dataset({'img': (('date', 'row', 'col'), (truth + aurora)[None]), 'sza': (('date', 'row', 'col'), sza[None]), 'dza': (('date', 'row', 'col'), np.full_like(x, 45)[None]), 'id': xr.DataArray('WIC')}, coords={'date': [np.datetime64('2001-01-01')], 'row': np.arange(20), 'col': np.arange(40)})
    return image, truth


SETTINGS = {'sigma_model_coefficients': [0, 400], 'sKnots': [-3.5, 0, 3.5], 'sOrder': 1, 'monotonic': True, 'stop': 1e-8}


def test_zero_strength_preserves_unconstrained_fit_exactly():
    image = _synthetic_wic().isel(date=[0])
    ordinary = background.backgroundmodel_BS(image, sigma_model_coefficients=[0, 100])
    disabled = background.backgroundmodel_BS(image, sigma_model_coefficients=[0, 100], background_prior=np.zeros_like(image.img), background_prior_strength=0)
    for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma']:
        np.testing.assert_array_equal(ordinary[name], disabled[name])
    assert 'background_prior_strength' not in ordinary.dgmodel.attrs


def test_all_auroral_frame_retains_prior_even_with_numerical_damping():
    image, truth = closeup(aurora=400)
    free = background.backgroundmodel_BS(image, **SETTINGS)
    anchored = background.backgroundmodel_BS(image, background_prior=truth[None], background_prior_strength=1, dampingVal=.01, **SETTINGS)
    assert np.mean(free.dgmodel.values[0] - truth) > 399
    np.testing.assert_allclose(anchored.dgmodel.values[0], truth, atol=1e-8)
    assert anchored.dgfit_converged.item()
    assert np.count_nonzero(anchored.dgweight) == 0
    assert anchored.dgmodel.attrs['background_prior_trace_fraction'] == 1
    assert np.nanmedian(anchored.dgimg) > 399


def test_quiet_frame_can_depart_from_prior_and_strength_controls_departure():
    image, truth = closeup(offset=20)
    prior = (truth - 20)[None]
    weak = background.backgroundmodel_BS(image, background_prior=prior, background_prior_strength=.1, **SETTINGS)
    strong = background.backgroundmodel_BS(image, background_prior=prior, background_prior_strength=10, **SETTINGS)
    weak_change = np.mean(weak.dgmodel.values - prior)
    strong_change = np.mean(strong.dgmodel.values - prior)
    assert 17 < weak_change < 20
    assert 0 < strong_change < 3
    assert weak.dgfit_converged.item() and strong.dgfit_converged.item()


def test_narrow_cubic_monotonic_prior_is_retained_when_all_data_are_rejected():
    image, truth = closeup(aurora=400)
    result = background.backgroundmodel_BS(image, background_prior=truth[None], background_prior_strength=1, dampingVal=.01, sigma_model_coefficients=[0, 400], monotonic=True, stop=1e-8)
    np.testing.assert_allclose(result.dgmodel.values[0], truth, atol=1e-6)
    assert np.count_nonzero(result.dgweight) == 0


def test_moderate_widespread_aurora_exposes_strength_tradeoff():
    image, truth = closeup(aurora=40)
    weak = background.backgroundmodel_BS(image, background_prior=truth[None], background_prior_strength=.1, **SETTINGS)
    strong = background.backgroundmodel_BS(image, background_prior=truth[None], background_prior_strength=10, **SETTINGS)
    assert np.median(weak.dgimg) < 10
    assert np.median(strong.dgimg) > 35
    # This is empirical shrinkage, not proof that the auroral signal is intact.
    assert np.median(strong.dgimg) < 40


def test_prior_normal_equations_match_direct_weighted_prediction_penalty():
    image, truth = closeup(offset=20)
    prior = (truth - 20)[None]
    image['img_variance'] = xr.ones_like(image.img)*np.linspace(1, 2, 40)
    result = background.backgroundmodel_BS(image, background_prior=prior, background_prior_strength=.7, variance_weighted=True, sigma_model_coefficients=[0, 1e20], sKnots=[-3.5, 0, 3.5], sOrder=1, stop=1e-10)
    # Both truth and prior are linear functions in the same spline space.
    np.testing.assert_allclose(result.dgmodel, (image.img + .7*prior)/1.7, atol=1e-8)


def test_stored_directional_amplitudes_reconstruct_default_model_exactly():
    image = _directional_background(True)
    options = {'spatial_model': 'legacy_x_azimuth', 'spencer_taper': True, 'sigma_model_coefficients': [0, 100]}
    ordinary = background.backgroundmodel_BS(image, **options)
    stored = background.backgroundmodel_BS(image, store_spatial_coefficients=True, **options)
    for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma']:
        np.testing.assert_array_equal(ordinary[name], stored[name])
    coefficients = stored.dgmodel_spatial_coefficients.values
    n_base = stored.dgmodel_spatial_coefficients.attrs['base_count']
    x = np.cos(np.deg2rad(image.sza.values))/np.cos(np.deg2rad(image.dza.values))
    knots = np.r_[np.repeat(stored.dgmodel.attrs['spatial_knots'][0], 3), stored.dgmodel.attrs['spatial_knots'], np.repeat(stored.dgmodel.attrs['spatial_knots'][-1], 3)]
    design = background.BSpline(knots, np.eye(n_base), 3)(x)
    prediction = np.einsum('trcs,ts->trc', design, coefficients[:, :n_base])
    prediction += np.einsum('trcs,ts->trc', design, coefficients[:, n_base:])*image.solar_view_sine_factor.values*background._spencer_taper_weights(image.sza.values)
    support = np.isfinite(stored.dgmodel.values)
    np.testing.assert_allclose(prediction[support], stored.dgmodel.values[support], atol=1e-10)


@pytest.mark.parametrize('strength', [-1, np.nan, np.inf, True])
def test_invalid_prior_strength_is_rejected(strength):
    image, truth = closeup()
    with pytest.raises(ValueError, match='background_prior_strength'):
        background.backgroundmodel_BS(image, background_prior=truth[None], background_prior_strength=strength, **SETTINGS)


def test_prior_shape_support_and_scope_are_validated():
    image, truth = closeup()
    for prior in [truth, np.full_like(image.img, np.nan)]:
        with pytest.raises(ValueError, match='background_prior must match'):
            background.backgroundmodel_BS(image, background_prior=prior, background_prior_strength=1, **SETTINGS)
    with pytest.raises(ValueError, match='supplied sigma_model_coefficients'):
        background.backgroundmodel_BS(image, background_prior=truth[None], background_prior_strength=1)
    with pytest.raises(ValueError, match='single frame'):
        background.backgroundmodel_BS(_synthetic_wic(), background_prior=np.zeros_like(_synthetic_wic().img), background_prior_strength=1, **SETTINGS)
    shifted = image.img.assign_coords(date=[np.datetime64('2001-01-02')])
    with pytest.raises(ValueError):
        background.backgroundmodel_BS(image, background_prior=shifted, background_prior_strength=1, **SETTINGS)


@pytest.mark.parametrize('fallback', [False, True])
def test_prior_projection_failure_honors_monotonic_fallback(monkeypatch, fallback):
    def failed_projection(*args, **kwargs):
        raise RuntimeError('Injected projection failure')

    image, truth = closeup(offset=20)
    monkeypatch.setattr(background, '_solve_monotonic_bs', failed_projection)
    with pytest.warns(RuntimeWarning, match='Injected projection failure'):
        result = background.backgroundmodel_BS(image, background_prior=(truth - 20)[None], background_prior_strength=1, monotonic_fallback=fallback, **SETTINGS)
    if fallback:
        assert result.dgsolver_flag.item() == 2
        assert result.dgfit_converged.item()
        assert np.all(np.isfinite(result.dgmodel))
    else:
        assert result.dgsolver_flag.item() == 3
        assert not result.dgfit_converged.item()
        assert np.all(np.isnan(result.dgmodel))
