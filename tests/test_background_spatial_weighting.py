import numpy as np
import pytest
import xarray as xr

from fuvpy.src.background import _bin_population_weights, backgroundmodel_BS


SPATIAL_BINS = np.arange(-0.2, 0.3 + 0.1, 0.1)


def test_pixels_share_bin_weight_and_bins_have_equal_influence():
    fraction = np.array([[-0.15, -0.14, -0.13, 0.05, 0.06, 0.25]])
    valid = np.ones_like(fraction, dtype=bool)

    weights = _bin_population_weights(fraction, valid, SPATIAL_BINS)

    np.testing.assert_allclose(weights[0, :3], weights[0, 0])
    np.testing.assert_allclose(weights[0, 3:5], weights[0, 3])

    influence = [
        np.sum(weights[0, :3]**2),
        np.sum(weights[0, 3:5]**2),
        weights[0, 5]**2,
    ]
    np.testing.assert_allclose(influence, influence[0])
    np.testing.assert_allclose(np.mean(weights[valid]**2), 1)


def test_bin_population_weights_are_calculated_per_frame():
    fraction = np.array([
        [0.01, 0.02, 0.11],
        [0.01, 0.11, 0.12],
    ])
    valid = np.ones_like(fraction, dtype=bool)

    weights = _bin_population_weights(fraction, valid, SPATIAL_BINS)

    assert weights[0, 0] == weights[0, 1]
    assert weights[0, 0] < weights[0, 2]
    assert weights[1, 1] == weights[1, 2]
    assert weights[1, 0] > weights[1, 1]
    np.testing.assert_allclose(np.mean(weights**2, axis=1), 1)


def test_empty_bins_and_out_of_fit_pixels_do_not_affect_weights():
    fraction = np.array([[0.01, 0.02, 0.11, 0.25, -0.15]])
    valid = np.array([[True, True, True, False, False]])

    weights = _bin_population_weights(fraction, valid, SPATIAL_BINS)
    reference = _bin_population_weights(
        fraction[:, :3],
        valid[:, :3],
        SPATIAL_BINS,
    )

    np.testing.assert_allclose(weights[0, :3], reference[0])
    assert np.all(np.isnan(weights[~valid]))


def test_bin_population_option_is_limited_to_legacy_fuv_coordinate():
    vis = xr.Dataset(
        data_vars={
            "id": "VIS",
            "sza": (("row", "col"), [[60.0]]),
            "img": (("row", "col"), [[1.0]]),
        },
    )

    with pytest.raises(ValueError, match="legacy FUV coordinate"):
        backgroundmodel_BS(vis, spatial_weighting="bin_population")
