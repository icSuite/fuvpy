import numpy as np

from fuvpy.src.background import backgroundmodel_BS
from test_background_dgsigma_freeze import _synthetic_wic


def test_fit_diagnostics_are_absent_by_default():
    result = backgroundmodel_BS(_synthetic_wic(), stop=1e-5)

    assert "dgfit_iterations" not in result
    assert "dgfit_converged" not in result
    assert result["dgmodel"].attrs["fit_iterations"] >= 2
    assert result["dgmodel"].attrs["fit_converged"] == 1


def test_fit_diagnostics_report_scalar_convergence():
    result = backgroundmodel_BS(_synthetic_wic(), stop=1e-5, diagnostics=True)

    assert result["dgfit_iterations"].dims == ()
    assert result["dgfit_converged"].dims == ()
    assert int(result["dgfit_iterations"]) >= 2
    assert bool(result["dgfit_converged"])
    assert np.issubdtype(result["dgfit_iterations"].dtype, np.integer)


def test_fit_configuration_is_stored_on_background_model():
    result = backgroundmodel_BS(_synthetic_wic(), sKnots=[-3.5, 0, 3.5], sOrder=2, n_tKnots=3, tOrder=1, dampingVal=1e-3, tukeyVal=4, noise_model="binned", spatial_weighting="bin_population", dgsigma_freeze_iteration=0, input_variance="unused_variance", variance_weighted=False, monotonic=True, monotonic_solver="normal", monotonic_fallback=True, stop=1e-5)
    attrs = result["dgmodel"].attrs

    np.testing.assert_array_equal(attrs["spatial_knots"], [-3.5, 0, 3.5])
    np.testing.assert_array_equal(attrs["temporal_knots_minutes_from_first_image"], [0, 2, 4])
    assert attrs["spatial_order"] == 2
    assert attrs["temporal_knot_count"] == 3
    assert attrs["temporal_order"] == 1
    assert attrs["noise_model"] == "binned"
    assert attrs["spatial_weighting"] == "bin_population"
    assert attrs["dgsigma_freeze_iteration"] == 0
    assert attrs["robust_weighting"] == "tukey_bisquare"
    assert attrs["tukey_value"] == 4
    assert attrs["input_variance"] == "unused_variance"
    assert attrs["variance_weighted"] == "false"
    assert attrs["monotonic_requested"] == "true"
    assert attrs["monotonic_solver"] == "normal"
    assert attrs["monotonic_fallback"] == "true"
    assert attrs["spencer_taper"] == "disabled"
    assert attrs["damping_value"] == 1e-3
    assert attrs["damping_scale"] > 0
    assert attrs["convergence_tolerance"] == 1e-5
    assert attrs["maximum_iterations"] == 100


def test_background_model_rejects_alternative_input_field():
    with np.testing.assert_raises_regex(TypeError, "inImg"):
        backgroundmodel_BS(_synthetic_wic(), inImg="other_image")
