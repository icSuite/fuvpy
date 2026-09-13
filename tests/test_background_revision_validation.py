import importlib.util
from pathlib import Path

import numpy as np
import pytest
import xarray as xr


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "validate_background_revision.py"
SPEC = importlib.util.spec_from_file_location("validate_background_revision", SCRIPT_PATH)
validation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validation)


def test_variants_are_the_six_cumulative_step_7_comparisons():
    names = [name for name, settings in validation.VARIANTS]
    assert names == [
        "released",
        "corrected_dgsigma",
        "bin_population",
        "frozen_dgsigma",
        "measurement_weighted",
        "monotonic"
    ]

    settings = dict(validation.VARIANTS)
    assert settings["released"] == {}
    assert settings["corrected_dgsigma"] == {"noise_model": "binned"}
    assert settings["bin_population"]["spatial_weighting"] == "bin_population"
    assert settings["frozen_dgsigma"]["dgsigma_freeze_iteration"] == 0
    assert settings["measurement_weighted"]["variance_weighted"] is True
    assert settings["monotonic"]["monotonic_solver"] == "normal"
    assert settings["monotonic"]["monotonic_fallback"] is False

    previous = {}
    for name, current in validation.VARIANTS:
        for key, value in previous.items():
            assert current[key] == value
        previous = current


def test_representative_frames_include_orbit_endpoints():
    np.testing.assert_array_equal(validation._representative_frames(10, 3), [0, 4, 9])
    np.testing.assert_array_equal(validation._representative_frames(2, 5), [0, 1])


def test_plot_limits_are_shared_and_robust():
    results = {
        "first": xr.Dataset({"dgimg": (("date", "row", "col"), [[[-2.0, 1.0]]])}),
        "second": xr.Dataset({"dgimg": (("date", "row", "col"), [[[2.0, np.nan]]])})
    }

    lower, upper = validation._plot_limits(results, "dgimg", frame=0, symmetric=True)

    assert lower == -upper
    assert lower < 0 < upper


def _variance_dataset(sensor="WIC", include_variance=False):
    dataset = xr.Dataset(
        {
            "img": (
                ("date", "row", "col"),
                np.array([[[0.0, 10.0]], [[20.0, np.nan]]])
            ),
            "sza": (
                ("date", "row", "col"),
                np.full((2, 1, 2), 60.0)
            ),
            "dza": (
                ("date", "row", "col"),
                np.zeros((2, 1, 2))
            ),
            "id": xr.DataArray(sensor)
        },
        coords={
            "date": [np.datetime64("2000-01-01"), np.datetime64("2000-01-02")],
            "row": [0],
            "col": [0, 1]
        }
    )
    if include_variance:
        dataset["img_variance"] = (
            dataset["img"].dims,
            np.array([[[7.0, 8.0]], [[9.0, np.nan]]])
        )
        dataset["img_variance"].attrs["poisson_offset"] = 0.5
        if sensor == "WIC":
            dataset["img_variance"].attrs["calculated_from"] = "img before optional WIC reflattening"

    return dataset


def test_stored_wic_variance_is_preserved():
    dataset = _variance_dataset(include_variance=True)
    original = dataset["img_variance"].copy(deep=True)

    result = validation._attach_measurement_variance(dataset, "WIC")

    np.testing.assert_allclose(result["img_variance"], original)
    assert result["img_variance"].attrs == original.attrs


def test_wic_without_stored_pre_reflattening_variance_is_rejected():
    with pytest.raises(ValueError, match="read_idl"):
        validation._attach_measurement_variance(_variance_dataset(), "WIC")


def test_invalid_stored_variance_on_fitted_support_is_rejected():
    dataset = _variance_dataset(include_variance=True)
    dataset["img_variance"].values[0, 0, 0] = 0.0

    with pytest.raises(ValueError, match="strictly positive"):
        validation._attach_measurement_variance(dataset, "WIC")


@pytest.mark.parametrize("poisson_offset", [None, 0.0, np.nan])
def test_stored_variance_requires_half_count_offset(poisson_offset):
    dataset = _variance_dataset(include_variance=True)
    if poisson_offset is None:
        del dataset["img_variance"].attrs["poisson_offset"]
    else:
        dataset["img_variance"].attrs["poisson_offset"] = poisson_offset

    with pytest.raises(ValueError, match="poisson_offset=0.5"):
        validation._attach_measurement_variance(dataset, "WIC")


@pytest.mark.parametrize(
    "calculated_from",
    [None, "img after optional WIC reflattening"]
)
def test_stored_wic_variance_requires_pre_reflattening_provenance(calculated_from):
    dataset = _variance_dataset(include_variance=True)
    if calculated_from is None:
        del dataset["img_variance"].attrs["calculated_from"]
    else:
        dataset["img_variance"].attrs["calculated_from"] = calculated_from

    with pytest.raises(ValueError, match="before optional WIC reflattening"):
        validation._attach_measurement_variance(dataset, "WIC")


def test_missing_si_variance_is_attached_from_calibration(monkeypatch):
    dataset = _variance_dataset(sensor="SI12")
    multiplier = np.full(dataset["img"].shape, 2.0)
    monkeypatch.setattr(validation, "calibration_multiplier", lambda dates, sensor: multiplier)

    result = validation._attach_measurement_variance(dataset, "SI12")
    expected = 2*dataset["img"] + 0.5*2**2

    assert result["img_variance"].dims == dataset["img"].dims
    np.testing.assert_allclose(result["img_variance"], expected)
    assert result["img_variance"].attrs["poisson_offset"] == 0.5


def test_cli_parser_requires_explicit_sensor_paths_and_output():
    arguments = validation.build_parser().parse_args(["--wic", "wic.nc", "--si12", "si12.nc", "--si13", "si13.nc", "--output-dir", "output"])

    assert arguments.wic == Path("wic.nc")
    assert arguments.si12 == Path("si12.nc")
    assert arguments.si13 == Path("si13.nc")
    assert arguments.output_dir == Path("output")
