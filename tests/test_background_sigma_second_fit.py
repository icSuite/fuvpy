import numpy as np
import pytest

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic


def test_second_estimate_uses_robust_weights_and_then_freezes(monkeypatch):
    calls = []
    real_estimator = background._noisemodel_binned

    def record(fraction, image, model, weights, knots, **kwargs):
        result = real_estimator(fraction, image, model, weights, knots, **kwargs)
        calls.append({'weights': weights.copy(), 'sigma': result[0].copy(), 'coefficients': result[1].copy()})
        return result

    monkeypatch.setattr(background, '_noisemodel_binned', record)
    result = background.backgroundmodel_BS(_synthetic_wic(), store_sigma_model=True, dgsigma_freeze_iteration=2, stop=1e-5)
    assert len(calls) == 2
    assert np.any(calls[1]['weights'] < calls[0]['weights'])
    np.testing.assert_array_equal(result.dgsigma.values.ravel(), calls[1]['sigma'])
    np.testing.assert_array_equal(result.dgmodel.attrs['sigma_model_coefficients'], calls[1]['coefficients'])
    expected = (1 - np.minimum(1, (result.dgimg.values/(5*result.dgsigma.values))**2))**2
    np.testing.assert_allclose(result.dgweight.values, expected)
    assert result.dgfit_iterations.item() >= 3
    assert result.dgmodel.attrs['dgsigma_estimation'] == 'first_two_fits_then_frozen'


def test_final_estimate_is_consumed_even_with_loose_stopping_tolerance(monkeypatch):
    estimates = []

    def changing(fraction, image, model, weights, knots):
        estimates.append(100. if not estimates else .5)
        return np.full_like(image, estimates[-1])

    monkeypatch.setattr(background, '_noisemodel_binned', changing)
    result = background.backgroundmodel_BS(_synthetic_wic(), dgsigma_freeze_iteration=2, stop=1e5)
    assert estimates == [100., .5]
    assert result.dgfit_iterations.item() == 3
    np.testing.assert_array_equal(result.dgsigma.values, .5)


def test_explicit_default_is_bitwise_identical():
    implicit = background.backgroundmodel_BS(_synthetic_wic(), stop=1e-5)
    explicit = background.backgroundmodel_BS(_synthetic_wic(), dgsigma_freeze_iteration=1, stop=1e-5)
    for name in ['dgmodel', 'dgimg', 'dgsigma', 'dgweight']:
        np.testing.assert_array_equal(implicit[name].values, explicit[name].values)
    np.testing.assert_equal(implicit.dgmodel.attrs, explicit.dgmodel.attrs)


@pytest.mark.parametrize('value', [True, 0, 3, 1.5])
def test_invalid_freeze_iteration(value):
    with pytest.raises(ValueError, match='dgsigma_freeze_iteration'):
        background.backgroundmodel_BS(_synthetic_wic(), dgsigma_freeze_iteration=value)


def test_supplied_relation_cannot_be_reestimated():
    with pytest.raises(ValueError, match='supplied sigma_model_coefficients'):
        background.backgroundmodel_BS(_synthetic_wic(), sigma_model_coefficients=[0, 4], dgsigma_freeze_iteration=2)
