import numpy as np
import pytest

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic


def _with_variance(imgs, value=1.0):
    imgs = imgs.copy(deep=True)
    imgs["img_variance"] = (
        imgs["img"].dims,
        np.full_like(imgs["img"].values, value),
    )
    return imgs


def test_constant_variance_matches_unweighted_fit_without_damping():
    imgs = _with_variance(_synthetic_wic(), value=4.0)

    ordinary = background.backgroundmodel_BS(
        imgs,
        dampingVal=0,
        stop=1e-5,
    )
    weighted = background.backgroundmodel_BS(
        imgs,
        dampingVal=0,
        stop=1e-5,
        input_variance="img_variance",
        variance_weighted=True,
    )

    for field in ["dgmodel", "dgimg", "dgweight", "dgsigma"]:
        np.testing.assert_allclose(
            weighted[field],
            ordinary[field],
            rtol=1e-10,
            atol=1e-10,
        )


def test_high_variance_observation_has_less_influence():
    clean = _synthetic_wic()
    contaminated = _with_variance(clean)
    target = dict(date=1, row=2, col=3)
    contaminated["img"].values[1, 2, 3] += 100.0
    contaminated["img_variance"].values[1, 2, 3] = 1e6

    baseline = background.backgroundmodel_BS(
        clean,
        dampingVal=0,
        tukeyVal=1e6,
        stop=1e-5,
    )
    ordinary = background.backgroundmodel_BS(
        contaminated,
        dampingVal=0,
        tukeyVal=1e6,
        stop=1e-5,
    )
    weighted = background.backgroundmodel_BS(
        contaminated,
        dampingVal=0,
        tukeyVal=1e6,
        stop=1e-5,
        input_variance="img_variance",
        variance_weighted=True,
    )

    baseline_value = float(baseline["dgmodel"].isel(**target))
    ordinary_error = abs(
        float(ordinary["dgmodel"].isel(**target)) - baseline_value
    )
    weighted_error = abs(
        float(weighted["dgmodel"].isel(**target)) - baseline_value
    )

    assert weighted_error < ordinary_error


@pytest.mark.parametrize(
    "noise_model, function_name",
    [
        ("legacy", "_noisemodel"),
        ("binned", "_noisemodel_binned"),
    ],
)
def test_measurement_variance_is_not_passed_to_dgsigma(
    monkeypatch,
    noise_model,
    function_name,
):
    captured_noise_weights = []

    def capture_noisemodel(fraction, image, model, weights, knots):
        captured_noise_weights.append(weights.copy())
        return np.full_like(image, 1000.0)

    monkeypatch.setattr(background, function_name, capture_noisemodel)
    imgs = _with_variance(_synthetic_wic())
    imgs["img_variance"].values[:, :, :4] = 100.0

    background.backgroundmodel_BS(
        imgs,
        dgsigma_freeze_iteration=0,
        noise_model=noise_model,
        stop=1e-5,
    )
    background.backgroundmodel_BS(
        imgs,
        dgsigma_freeze_iteration=0,
        noise_model=noise_model,
        stop=1e-5,
        input_variance="img_variance",
        variance_weighted=True,
    )

    assert len(captured_noise_weights) == 2
    np.testing.assert_array_equal(
        captured_noise_weights[1],
        captured_noise_weights[0],
    )


@pytest.mark.parametrize(
    "value, message",
    [
        (np.nan, "finite"),
        (np.inf, "finite"),
        (0.0, "strictly positive"),
        (-1.0, "strictly positive"),
    ],
)
def test_invalid_variance_at_fitted_pixel_is_rejected(value, message):
    imgs = _with_variance(_synthetic_wic())
    imgs["img_variance"].values[0, 0, 0] = value

    with pytest.raises(ValueError, match=message):
        background.backgroundmodel_BS(
            imgs,
            input_variance="img_variance",
            variance_weighted=True,
        )


def test_variance_weighting_requires_named_dataset_variable():
    with pytest.raises(ValueError, match="must name a Dataset variable"):
        background.backgroundmodel_BS(
            _synthetic_wic(),
            variance_weighted=True,
        )
