import numpy as np
import pytest
import xarray as xr

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic
from test_background_measurement_weighting import _with_variance


def test_asymmetric_coordinate_preserves_dayside_and_compresses_nightside():
    legacy = np.array([-2.0, -0.2, -0.05, 0.0, 0.4, 2.0])
    width = 0.2

    coordinate = background._asymmetric_legacy_coordinate(legacy, width)
    expected = np.where(legacy >= 0, legacy, width*np.tanh(legacy/width))

    np.testing.assert_allclose(coordinate, expected)
    np.testing.assert_array_equal(coordinate[legacy >= 0], legacy[legacy >= 0])


def test_default_legacy_fit_is_unchanged_when_asymmetry_is_disabled():
    imgs = _synthetic_wic()

    default = background.backgroundmodel_BS(imgs, stop=1e-5)
    explicit = background.backgroundmodel_BS(imgs, asymmetric_legacy=False, stop=1e-5)

    for field in ["dgmodel", "dgimg", "dgweight", "dgsigma"]:
        np.testing.assert_array_equal(explicit[field], default[field])


def test_default_width_and_knots_match_the_experimental_coordinate():
    width = 0.2
    conceptual_knots = np.array([-width, -0.75*width, -0.25*width, 0, 0.25, 1.5, 3.5])
    expected = conceptual_knots.copy()
    padding = 1e-6*(expected[-1] - expected[0])
    expected[0] -= padding
    expected[-1] += padding

    np.testing.assert_allclose(background._asymmetric_legacy_knots(width), expected)

    imgs = _synthetic_wic()
    default_width = background.backgroundmodel_BS(imgs, asymmetric_legacy=True, stop=1e-5)
    explicit_width = background.backgroundmodel_BS(imgs, asymmetric_legacy=True, nightside_width=0.2, stop=1e-5)
    np.testing.assert_array_equal(default_width["dgmodel"], explicit_width["dgmodel"])


def test_custom_nightside_width_changes_only_negative_coordinate_values():
    legacy = np.array([-1.0, -0.1, 0.0, 0.25, 1.0])
    coordinate = background._asymmetric_legacy_coordinate(legacy, 0.4)

    np.testing.assert_allclose(coordinate[:2], 0.4*np.tanh(legacy[:2]/0.4))
    np.testing.assert_array_equal(coordinate[2:], legacy[2:])


@pytest.mark.parametrize("width", [0, -0.1, np.nan, np.inf, None, "invalid"])
def test_nightside_width_must_be_finite_and_positive(width):
    with pytest.raises(ValueError, match="finite number greater than zero"):
        background.backgroundmodel_BS(xr.Dataset(), nightside_width=width)


def test_asymmetric_coordinate_is_rejected_for_vis():
    vis = xr.Dataset({"id": "VIS", "sza": (("row", "col"), [[60.0]]), "img": (("row", "col"), [[1.0]])})

    with pytest.raises(ValueError, match="not available for VIS"):
        background.backgroundmodel_BS(vis, asymmetric_legacy=True)


def test_custom_spatial_knots_are_respected(monkeypatch):
    def unexpected_default_knots(width):
        raise AssertionError("default asymmetric knots should not be requested")

    monkeypatch.setattr(background, "_asymmetric_legacy_knots", unexpected_default_knots)
    custom_knots = [-0.3, -0.15, -0.05, 0, 0.25, 1.5, 3.5]

    result = background.backgroundmodel_BS(_synthetic_wic(), asymmetric_legacy=True, sKnots=custom_knots, stop=1e-5)
    assert np.any(np.isfinite(result["dgmodel"]))


def test_asymmetric_coordinate_combines_with_revised_fit_options():
    imgs = _with_variance(_synthetic_wic())
    imgs["img_variance"] = 1 + np.abs(imgs["img"])

    result = background.backgroundmodel_BS(imgs, asymmetric_legacy=True, noise_model="binned", spatial_weighting="bin_population", dgsigma_freeze_iteration=0, input_variance="img_variance", variance_weighted=True, monotonic=True, stop=1e-5)

    assert int(result["dgsolver_flag"]) == 1
    assert np.any(np.isfinite(result["dgmodel"]))
    assert np.any(np.isfinite(result["dgsigma"]))
