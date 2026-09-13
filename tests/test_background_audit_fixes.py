import numpy as np
import pytest

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic


def test_normal_matrix_scale_ignores_unsupported_columns():
    normal_matrix = np.diag([0.0, 1e-30, 2.0, 4.0])

    assert background._normal_matrix_scale(normal_matrix) == 3.0


def test_normal_matrix_scale_rejects_no_supported_columns():
    with pytest.raises(ValueError, match="no supported coefficients"):
        background._normal_matrix_scale(np.zeros((3, 3)))


def test_unknown_keyword_is_rejected():
    with pytest.raises(TypeError, match="nonsense_option"):
        background.backgroundmodel_BS(_synthetic_wic(), nonsense_option=True)


@pytest.mark.parametrize("option", ["dzalim", "minlat"])
def test_obsolete_geometry_keywords_are_rejected(option):
    with pytest.raises(TypeError, match=option):
        background.backgroundmodel_BS(_synthetic_wic(), **{option: 75})


def test_returned_fields_are_masked_outside_fit_support():
    imgs = _synthetic_wic()
    imgs["dza"].values[0, 0, 0] = 77.0
    imgs["img"].values[0, 0, 1] = np.nan

    result = background.backgroundmodel_BS(imgs, stop=1e-5)

    for field in ["dgmodel", "dgimg", "dgweight", "dgsigma"]:
        assert np.isnan(result[field].values[0, 0, 0])
        assert np.isnan(result[field].values[0, 0, 1])


def test_dgsigma_uses_tukey_times_squared_spatial_multiplier(monkeypatch):
    captured_weights = []

    def capture_noisemodel(fraction, image, model, weights, knots):
        captured_weights.append(weights.copy())
        return np.full_like(image, 1000.0)

    monkeypatch.setattr(background, "_noisemodel_binned", capture_noisemodel)
    imgs = _synthetic_wic()

    background.backgroundmodel_BS(imgs, stop=1e-5)

    fraction = np.cos(np.deg2rad(imgs["sza"].stack(z=("row", "col")).values))/np.cos(np.deg2rad(imgs["dza"].stack(z=("row", "col")).values))
    valid = (imgs["sza"].stack(z=("row", "col")).values >= 0) & np.isfinite(imgs["img"].stack(z=("row", "col")).values) & (fraction > -3.5) & (fraction < 3.5)
    spatial_bins = np.arange(-3.5, 3.5 + 0.1, 0.1)
    spatial_multiplier = background._bin_population_weights(fraction, valid, spatial_bins)

    assert len(captured_weights) == 1
    np.testing.assert_allclose(captured_weights[0], spatial_multiplier[valid]**2)
