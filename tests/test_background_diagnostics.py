import numpy as np

from fuvpy.src.background import backgroundmodel_BS
from test_background_dgsigma_freeze import _synthetic_wic


def test_fit_returns_a_copy_with_diagnostics():
    original = _synthetic_wic()
    result = backgroundmodel_BS(original, stop=1e-5)

    assert result is not original
    assert "dgmodel" not in original
    assert result["dgfit_iterations"].dims == ()
    assert result["dgfit_converged"].dims == ()
    assert int(result["dgfit_iterations"]) >= 2
    assert bool(result["dgfit_converged"])
    assert np.issubdtype(result["dgfit_iterations"].dtype, np.integer)
    assert result["dgmodel"].attrs["fit_iterations"] >= 2
    assert result["dgmodel"].attrs["fit_converged"] == 1


def test_fit_configuration_is_stored_on_background_model():
    result = backgroundmodel_BS(_synthetic_wic(), sKnots=[-3.5, 0, 3.5], sOrder=2, n_tKnots=3, tOrder=1, dampingVal=1e-3, tukeyVal=4, variance_weighted=False, monotonic=True, monotonic_fallback=True, stop=1e-5)
    attrs = result["dgmodel"].attrs

    np.testing.assert_array_equal(attrs["spatial_knots"], [-3.5, 0, 3.5])
    np.testing.assert_array_equal(attrs["temporal_knots_minutes_from_first_image"], [0, 2, 4])
    assert attrs["spatial_order"] == 2
    assert attrs["temporal_knot_count"] == 3
    assert attrs["temporal_order"] == 1
    assert attrs["noise_model"] == "binned"
    assert attrs["spatial_weighting"] == "bin_population"
    assert attrs["dgsigma_estimation"] == "initial_fit"
    assert attrs["robust_weighting"] == "tukey_bisquare"
    assert attrs["tukey_value"] == 4
    assert attrs["measurement_variance"] == "not_used"
    assert attrs["variance_weighted"] == "false"
    assert attrs["monotonic_requested"] == "true"
    assert attrs["monotonic_solver"] == "normal"
    assert attrs["monotonic_fallback"] == "true"
    assert attrs["spencer_taper"] == "disabled"
    assert attrs["damping_value"] == 1e-3
    assert attrs["damping_scale"] > 0
    assert attrs["convergence_tolerance"] == 1e-5
    assert attrs["maximum_iterations"] == 100
