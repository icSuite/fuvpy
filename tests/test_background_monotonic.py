from types import SimpleNamespace

import numpy as np
import pytest

import fuvpy.src.background as background
from test_background_dgsigma_freeze import _synthetic_wic
from test_background_measurement_weighting import _with_variance


def test_monotonic_transform_preserves_temporal_ordering():
    transform = background._monotonic_bs_transform(n_scp=3, n_tcp=2)
    free_parameters = np.array([10.0, -2.0, 1.0, 3.0, 2.0, 4.0])

    coefficients = (transform @ free_parameters).reshape((3, 2))

    np.testing.assert_allclose(coefficients, [[10.0, -2.0], [11.0, 1.0], [13.0, 5.0]])
    assert np.all(np.diff(coefficients, axis=0) >= 0)


def test_monotonic_fit_is_nondecreasing_in_legacy_coordinate():
    imgs = _synthetic_wic()
    result = background.backgroundmodel_BS(imgs, monotonic=True, stop=1e-5)

    fraction = (np.cos(np.deg2rad(imgs["sza"]))/np.cos(np.deg2rad(imgs["dza"]))).values

    np.testing.assert_array_equal(result["date"], imgs["date"])
    assert int(result["dgsolver_flag"]) == 1
    for frame in range(len(result["date"])):
        order = np.argsort(fraction[frame].ravel())
        ordered_model = result["dgmodel"].values[frame].ravel()[order]
        assert np.all(np.diff(ordered_model) >= -1e-10)


def test_normal_and_full_solvers_give_equivalent_weighted_models():
    imgs = _with_variance(_synthetic_wic())
    imgs["img_variance"] = 1 + np.abs(imgs["img"])
    options = {
        "monotonic": True,
        "spatial_weighting": "bin_population",
        "variance_weighted": True,
        "input_variance": "img_variance",
        "noise_model": "binned",
        "dgsigma_freeze_iteration": 0,
        "dampingVal": 1e-3,
        "stop": 1e-5
    }

    normal = background.backgroundmodel_BS(imgs, monotonic_solver="normal", **options)
    full = background.backgroundmodel_BS(imgs, monotonic_solver="full", **options)

    np.testing.assert_allclose(normal["dgmodel"], full["dgmodel"], rtol=1e-7, atol=1e-7)


def test_explicitly_disabling_monotonic_fit_matches_default():
    imgs = _synthetic_wic()

    ordinary = background.backgroundmodel_BS(imgs, stop=1e-5)
    explicit = background.backgroundmodel_BS(imgs, monotonic=False, monotonic_solver="full", stop=1e-5)

    for field in ["dgmodel", "dgimg", "dgweight", "dgsigma"]:
        np.testing.assert_array_equal(explicit[field], ordinary[field])
    assert int(ordinary["dgsolver_flag"]) == 0
    assert int(explicit["dgsolver_flag"]) == 0
    np.testing.assert_array_equal(ordinary["dgsolver_flag"].attrs["flag_values"], [0, 1, 2, 3, 4])
    assert ordinary["dgsolver_flag"].attrs["flag_meanings"] == "standard_requested monotonic_succeeded monotonic_failed_standard_fallback monotonic_failed_no_fallback monotonic_failed_standard_fallback_failed"


def test_monotonic_failure_warns_and_falls_back_once(monkeypatch):
    calls = 0

    def fail_once(*args, **kwargs):
        nonlocal calls
        calls += 1
        raise RuntimeError("forced monotonic failure")

    monkeypatch.setattr(background, "_solve_monotonic_bs", fail_once)

    with pytest.warns(RuntimeWarning, match="all remaining IRLS iterations"):
        result = background.backgroundmodel_BS(_synthetic_wic(), monotonic=True, monotonic_fallback=True, diagnostics=True, stop=1e-5)

    assert calls == 1
    assert int(result["dgsolver_flag"]) == 2
    assert result["dgmodel"].attrs["monotonic_constraint"] == "none"
    assert bool(result["dgfit_converged"])


def test_monotonic_failure_returns_flagged_nan_product_when_fallback_is_disabled(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("forced monotonic failure")

    monkeypatch.setattr(background, "_solve_monotonic_bs", fail)

    with pytest.warns(RuntimeWarning, match="all-NaN background product"):
        result = background.backgroundmodel_BS(_synthetic_wic(), monotonic=True, monotonic_fallback=False)

    assert int(result["dgsolver_flag"]) == 3
    assert not bool(result["dgfit_converged"])
    assert result["dgmodel"].attrs["monotonic_constraint"] == "not_available"
    assert "forced monotonic failure" in result["dgsolver_flag"].attrs["solver_message"]
    for field in ["dgmodel", "dgimg", "dgweight", "dgsigma"]:
        assert np.all(np.isnan(result[field]))


def test_failed_ordinary_fallback_returns_flagged_nan_product(monkeypatch):
    def fail_monotonic(*args, **kwargs):
        raise RuntimeError("forced monotonic failure")

    def fail_fallback(*args, **kwargs):
        raise np.linalg.LinAlgError("forced fallback failure")

    monkeypatch.setattr(background, "_solve_monotonic_bs", fail_monotonic)
    monkeypatch.setattr(background, "lstsq", fail_fallback)

    with pytest.warns(RuntimeWarning) as warnings_record:
        result = background.backgroundmodel_BS(_synthetic_wic(), monotonic=True, monotonic_fallback=True)

    assert len(warnings_record) == 2
    assert int(result["dgsolver_flag"]) == 4
    assert not bool(result["dgfit_converged"])
    assert np.all(np.isnan(result["dgmodel"]))
    assert "forced fallback failure" in result["dgsolver_flag"].attrs["solver_message"]


def test_monotonic_solver_failure_is_explicit(monkeypatch):
    failure = SimpleNamespace(
        success=False,
        optimality=1.0,
        message="forced failure",
        x=np.zeros(4)
    )
    monkeypatch.setattr(background, "lsq_linear", lambda *args, **kwargs: failure)

    with pytest.raises(RuntimeError, match="forced failure"):
        background._solve_monotonic_bs(weighted_design=np.eye(4), weighted_data=np.ones(4), damping=np.zeros((4, 4)), n_scp=2, n_tcp=2, solver="full")


def test_monotonic_solver_exception_becomes_runtime_error(monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("invalid bounded problem")

    monkeypatch.setattr(background, "lsq_linear", fail)

    with pytest.raises(RuntimeError, match="invalid bounded problem"):
        background._solve_monotonic_bs(weighted_design=np.eye(4), weighted_data=np.ones(4), damping=np.zeros((4, 4)), n_scp=2, n_tcp=2, solver="full")


def test_monotonic_solver_name_is_validated():
    with pytest.raises(ValueError, match="monotonic_solver"):
        background.backgroundmodel_BS(_synthetic_wic(), monotonic=True, monotonic_solver="unknown")
