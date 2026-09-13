import numpy as np
from fuvpy.src.background import _bin_population_weights


SPATIAL_BINS = np.arange(-0.2, 0.3 + 0.1, 0.1)


def test_bin_population_weights_equalize_bins_independently_per_frame():
    fraction = np.array([[-0.15, -0.14, -0.13, 0.05, 0.06, 0.25]])
    valid = np.ones_like(fraction, dtype=bool)

    weights = _bin_population_weights(fraction, valid, SPATIAL_BINS)

    np.testing.assert_allclose(weights[0, :3], weights[0, 0])
    np.testing.assert_allclose(weights[0, 3:5], weights[0, 3])

    influence = [
        np.sum(weights[0, :3]**2),
        np.sum(weights[0, 3:5]**2),
        weights[0, 5]**2
    ]
    np.testing.assert_allclose(influence, influence[0])
    np.testing.assert_allclose(np.mean(weights[valid]**2), 1)

    fraction = np.array([[0.01, 0.02, 0.11], [0.01, 0.11, 0.12]])
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
    reference = _bin_population_weights(fraction[:, :3], valid[:, :3], SPATIAL_BINS)

    np.testing.assert_allclose(weights[0, :3], reference[0])
    assert np.all(np.isnan(weights[~valid]))
