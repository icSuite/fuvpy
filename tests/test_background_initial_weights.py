import numpy as np
import pytest
import xarray as xr

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic


def test_all_ones_initialization_preserves_default_fields_and_damping_exactly():
    image = _synthetic_wic().isel(date=[0])
    options = {'sigma_model_coefficients': [0, 100], 'dampingVal': .01, 'stop': 1e-6}
    ordinary = background.backgroundmodel_BS(image, **options)
    seeded = background.backgroundmodel_BS(image, initial_robust_weights=np.ones_like(image.img), **options)
    for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma']:
        np.testing.assert_array_equal(ordinary[name], seeded[name])
    assert ordinary.dgmodel.attrs['damping_scale'] == seeded.dgmodel.attrs['damping_scale']


def test_initial_auroral_weights_prevent_nightside_background_capture():
    x = np.broadcast_to(np.linspace(-.9, -.1, 40), (20, 40))
    truth = 200 + 40*x
    aurora = np.zeros_like(x)
    aurora[4:] = 150
    sza = np.rad2deg(np.arccos(x*np.cos(np.deg2rad(45))))
    image = xr.Dataset({'img': (('date', 'row', 'col'), (truth + aurora)[None]), 'sza': (('date', 'row', 'col'), sza[None]), 'dza': (('date', 'row', 'col'), np.full_like(x, 45)[None]), 'id': xr.DataArray('WIC')}, coords={'date': [np.datetime64('2001-01-01')], 'row': np.arange(20), 'col': np.arange(40)})
    initial = (aurora == 0).astype(float)[None]
    initial[0, 0] = 0  # Quiet pixels mistakenly rejected by the first-stage fit.
    options = {'sigma_model_coefficients': [0, 400], 'sKnots': [-3.5, 0, 3.5], 'sOrder': 1, 'monotonic': True, 'stop': 1e-8}
    free = background.backgroundmodel_BS(image, **options)
    seeded = background.backgroundmodel_BS(image, initial_robust_weights=initial, **options)
    assert np.mean(free.dgmodel.values[0] - truth) > 149
    np.testing.assert_allclose(seeded.dgmodel.values[0], truth, atol=1e-8)
    assert np.min(seeded.dgweight.values[0, 0]) > .99
    assert np.median(seeded.dgimg.values[0, 4:]) > 149


def test_seed_changes_only_first_solve_and_not_damping_scale(monkeypatch):
    calls = []
    original = background._accumulate_normal_system

    def recording(*args):
        calls.append(args[5].copy())
        return original(*args)

    image = _synthetic_wic().isel(date=[0])
    initial = np.full_like(image.img, .3)
    options = {'sigma_model_coefficients': [0, 100], 'dampingVal': .01, 'stop': 1e-6}
    ordinary = background.backgroundmodel_BS(image, **options)
    monkeypatch.setattr(background, '_accumulate_normal_system', recording)
    seeded = background.backgroundmodel_BS(image, initial_robust_weights=initial, **options)
    np.testing.assert_array_equal(calls[0], np.ones_like(calls[0]))  # Damping reference.
    np.testing.assert_array_equal(calls[1], initial.reshape(-1))
    assert np.any(calls[2] != initial.reshape(-1))
    assert seeded.dgmodel.attrs['damping_scale'] == ordinary.dgmodel.attrs['damping_scale']


@pytest.mark.parametrize('value', [-.1, 1.1, np.nan, 0])
def test_invalid_initial_weights_are_rejected(value):
    image = _synthetic_wic().isel(date=[0])
    with pytest.raises(ValueError, match='initial_robust_weights'):
        background.backgroundmodel_BS(image, initial_robust_weights=np.full_like(image.img, value), sigma_model_coefficients=[0, 100])


def test_initial_weight_shape_scope_and_coordinates_are_checked():
    image = _synthetic_wic().isel(date=[0])
    with pytest.raises(ValueError, match='match img shape'):
        background.backgroundmodel_BS(image, initial_robust_weights=np.ones((6, 8)), sigma_model_coefficients=[0, 100])
    with pytest.raises(ValueError, match='single frame'):
        background.backgroundmodel_BS(_synthetic_wic(), initial_robust_weights=np.ones_like(_synthetic_wic().img), sigma_model_coefficients=[0, 100])
    with pytest.raises(ValueError, match='requires supplied sigma'):
        background.backgroundmodel_BS(image, initial_robust_weights=np.ones_like(image.img))
    with pytest.raises(ValueError, match='cannot combine'):
        background.backgroundmodel_BS(image, initial_robust_weights=np.ones_like(image.img), background_prior=image.img, background_prior_strength=1, sigma_model_coefficients=[0, 100])
    with pytest.raises(ValueError):
        background.backgroundmodel_BS(image, initial_robust_weights=image.img.assign_coords(date=[np.datetime64('2002-01-01')]), sigma_model_coefficients=[0, 100])
