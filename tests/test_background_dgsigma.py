import numpy as np

from fuvpy.src.background import _binned_noise_statistics, _noise_bins, _noisemodel_binned


SPATIAL_KNOTS = np.array([-3.5, -0.25, 0, 0.25, 1.5, 3.5])


def test_noise_bins_have_one_nightside_and_fourteen_dayside_bins():
    bins = _noise_bins(SPATIAL_KNOTS)

    assert bins[0] == -3.5
    assert bins[1] == 0
    assert bins[-1] == 3.5
    assert len(bins) - 1 == 15
    np.testing.assert_allclose(np.diff(bins[1:]), 0.25)


def test_binned_noise_model_recovers_known_nonnegative_relation():
    slope = 0.4
    intercept = 2.5
    fraction = np.array([-0.5, 0.125, 0.375, 0.625, 0.875])
    model = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    residual = np.sqrt(slope*model**2 + intercept)
    image = model + residual
    weights = np.ones_like(model)

    sigma = _noisemodel_binned(fraction, image, model, weights, SPATIAL_KNOTS)

    np.testing.assert_allclose(sigma, np.sqrt(slope*model**2 + intercept), rtol=1e-7, atol=1e-7)


def test_binned_noise_fit_is_invariant_to_one_bin_population():
    fraction = np.array([-0.5, 0.125, 0.375, 0.625])
    model_squared = np.array([1.0, 4.0, 9.0, 16.0])
    residual_squared = np.array([3.0, 6.0, 10.0, 14.0])
    model = np.sqrt(model_squared)
    image = model + np.sqrt(residual_squared)
    weights = np.ones_like(model)

    reference = _noisemodel_binned(fraction, image, model, weights, SPATIAL_KNOTS)

    repetitions = 100
    expanded_fraction = np.r_[fraction, np.repeat(fraction[1], repetitions)]
    expanded_model = np.r_[model, np.repeat(model[1], repetitions)]
    expanded_image = np.r_[image, np.repeat(image[1], repetitions)]
    expanded_weights = np.ones_like(expanded_model)

    expanded = _noisemodel_binned(expanded_fraction, expanded_image, expanded_model, expanded_weights, SPATIAL_KNOTS)
    x_bin, y_bin = _binned_noise_statistics(expanded_fraction, expanded_image, expanded_model, expanded_weights, SPATIAL_KNOTS)

    assert x_bin.size == 4
    assert y_bin.size == 4
    np.testing.assert_allclose(expanded[:model.size], reference)
