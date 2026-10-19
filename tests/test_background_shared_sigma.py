import numpy as np
import pytest
import xarray as xr

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic


def test_storing_sigma_relation_preserves_default_fit_exactly():
    images = _synthetic_wic()
    ordinary = background.backgroundmodel_BS(images, stop=1e-5)
    stored = background.backgroundmodel_BS(images, stop=1e-5, store_sigma_model=True)
    for name in ['dgmodel', 'dgimg', 'dgweight', 'dgsigma']:
        np.testing.assert_array_equal(stored[name].values, ordinary[name].values)
    coefficients = stored.dgmodel.attrs['sigma_model_coefficients']
    assert coefficients.shape == (2,)
    assert np.all(coefficients >= 0)
    assert 'sigma_model_coefficients' not in ordinary.dgmodel.attrs


def test_supplied_relation_skips_bins_and_reinitializes_robust_weights(monkeypatch):
    def no_binned_estimate(*args, **kwargs):
        pytest.fail('A supplied sigma relation must not estimate frame-bin statistics')

    monkeypatch.setattr(background, '_noisemodel_binned', no_binned_estimate)
    image = _synthetic_wic().isel(date=[0])
    reference = background.backgroundmodel_BS(image, sigma_model_coefficients=[.002, 9], stop=1e-5)
    image['dgweight'] = xr.zeros_like(image.img)
    image['dgsigma'] = xr.full_like(image.img, 1e9)
    repeated = background.backgroundmodel_BS(image, sigma_model_coefficients=[.002, 9], stop=1e-5)
    for name in ['dgmodel', 'dgweight', 'dgsigma']:
        np.testing.assert_array_equal(repeated[name].values, reference[name].values)
    np.testing.assert_allclose(reference.dgsigma**2, .002*reference.dgmodel**2 + 9)
    assert reference.dgmodel.attrs['temporal_order'] == 0
    assert reference.dgmodel.attrs['noise_model'] == 'supplied_quadratic'


@pytest.mark.parametrize('coefficients', [[0, 0], [-1, 2], [1, -2], [np.nan, 1], [1, np.inf], [1], [1, 2, 3]])
def test_invalid_sigma_relations_are_rejected(coefficients):
    with pytest.raises(ValueError, match='sigma_model_coefficients'):
        background.backgroundmodel_BS(_synthetic_wic(), sigma_model_coefficients=coefficients)


def test_independent_frames_recover_abrupt_background_with_injected_aurora():
    generator = np.random.default_rng(9268)
    fraction = np.broadcast_to(np.linspace(-.8, 1.2, 40), (20, 40))
    dza = np.full_like(fraction, 45.)
    sza = np.rad2deg(np.arccos(fraction*np.cos(np.deg2rad(dza))))
    aurora = np.zeros_like(fraction)
    aurora[7:10, 10:20] = 150
    for frame, offset in enumerate([0, 70, -50]):
        truth = 200 + 80*fraction + offset
        image = xr.Dataset({'img': (('date', 'row', 'col'), (truth + aurora + generator.normal(0, 2, truth.shape))[None]), 'sza': (('date', 'row', 'col'), sza[None]), 'dza': (('date', 'row', 'col'), dza[None]), 'id': xr.DataArray('WIC')}, coords={'date': [np.datetime64('2001-01-01') + frame*np.timedelta64(2, 'm')], 'row': np.arange(20), 'col': np.arange(40)})
        result = background.backgroundmodel_BS(image, sigma_model_coefficients=[0, 4], sKnots=[-3.5, 0, 3.5], sOrder=1, monotonic=True, stop=1e-6)
        assert result.dgfit_converged.item()
        assert np.sqrt(np.mean((result.dgmodel.values[0] - truth)**2)) < 1
        assert np.median(result.dgimg.values[0][aurora > 0]) > 145


def test_sigma_coefficients_survive_netcdf_round_trip(tmp_path):
    fitted = background.backgroundmodel_BS(_synthetic_wic(), store_sigma_model=True)
    path = tmp_path/'sigma_model.nc'
    fitted.to_netcdf(path)
    with xr.open_dataset(path) as restored:
        np.testing.assert_array_equal(restored.dgmodel.attrs['sigma_model_coefficients'], fitted.dgmodel.attrs['sigma_model_coefficients'])
