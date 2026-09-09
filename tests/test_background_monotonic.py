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

    np.testing.assert_allclose(
        coefficients,
        [[10.0, -2.0], [11.0, 1.0], [13.0, 5.0]],
    )
    assert np.all(np.diff(coefficients, axis=0) >= 0)


def test_monotonic_fit_is_nondecreasing_in_legacy_coordinate():
    imgs = _synthetic_wic()
    result = background.backgroundmodel_BS(
        imgs,
        monotonic=True,
        stop=1e-5,
    )

    fraction = (
        np.cos(np.deg2rad(imgs["sza"]))
        / np.cos(np.deg2rad(imgs["dza"]))
    ).values

    np.testing.assert_array_equal(result["date"], imgs["date"])
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
        "stop": 1e-5,
    }

    normal = background.backgroundmodel_BS(
        imgs,
        monotonic_solver="normal",
        **options,
    )
    full = background.backgroundmodel_BS(
        imgs,
        monotonic_solver="full",
        **options,
    )

    np.testing.assert_allclose(
        normal["dgmodel"],
        full["dgmodel"],
        rtol=1e-7,
        atol=1e-7,
    )


def test_explicitly_disabling_monotonic_fit_matches_default():
    imgs = _synthetic_wic()

    ordinary = background.backgroundmodel_BS(imgs, stop=1e-5)
    explicit = background.backgroundmodel_BS(
        imgs,
        monotonic=False,
        monotonic_solver="full",
        stop=1e-5,
    )

    for field in ["dgmodel", "dgimg", "dgweight", "dgsigma"]:
        np.testing.assert_array_equal(explicit[field], ordinary[field])


def test_monotonic_solver_failure_is_explicit(monkeypatch):
    failure = SimpleNamespace(
        success=False,
        optimality=1.0,
        message="forced failure",
        x=np.zeros(4),
    )
    monkeypatch.setattr(background, "lsq_linear", lambda *args, **kwargs: failure)

    with pytest.raises(RuntimeError, match="forced failure"):
        background._solve_monotonic_bs(
            weighted_design=np.eye(4),
            weighted_data=np.ones(4),
            damping=np.zeros((4, 4)),
            n_scp=2,
            n_tcp=2,
            solver="full",
        )


def test_monotonic_solver_exception_becomes_runtime_error(monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("invalid bounded problem")

    monkeypatch.setattr(background, "lsq_linear", fail)

    with pytest.raises(RuntimeError, match="invalid bounded problem"):
        background._solve_monotonic_bs(
            weighted_design=np.eye(4),
            weighted_data=np.ones(4),
            damping=np.zeros((4, 4)),
            n_scp=2,
            n_tcp=2,
            solver="full",
        )


def test_monotonic_solver_name_is_validated():
    with pytest.raises(ValueError, match="monotonic_solver"):
        background.backgroundmodel_BS(
            _synthetic_wic(),
            monotonic=True,
            monotonic_solver="unknown",
        )
