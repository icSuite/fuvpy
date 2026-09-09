import importlib.util
from pathlib import Path

import numpy as np
import xarray as xr

from fuvpy.src.measurement_uncertainty import wic_variance_coefficient


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "scripts"
    / "validate_background_revision.py"
)
SPEC = importlib.util.spec_from_file_location(
    "validate_background_revision",
    SCRIPT_PATH,
)
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
        "monotonic",
    ]

    settings = dict(validation.VARIANTS)
    assert settings["released"] == {}
    assert settings["corrected_dgsigma"] == {"noise_model": "binned"}
    assert settings["bin_population"]["spatial_weighting"] == "bin_population"
    assert settings["frozen_dgsigma"]["dgsigma_freeze_iteration"] == 0
    assert settings["measurement_weighted"]["variance_weighted"] is True
    assert settings["monotonic"]["monotonic_solver"] == "normal"

    previous = {}
    for name, current in validation.VARIANTS:
        for key, value in previous.items():
            assert current[key] == value
        previous = current


def test_representative_frames_include_orbit_endpoints():
    np.testing.assert_array_equal(
        validation._representative_frames(10, 3),
        [0, 4, 9],
    )
    np.testing.assert_array_equal(
        validation._representative_frames(2, 5),
        [0, 1],
    )


def test_plot_limits_are_shared_and_robust():
    results = {
        "first": xr.Dataset(
            {"dgimg": (("date", "row", "col"), [[[-2.0, 1.0]]])}
        ),
        "second": xr.Dataset(
            {"dgimg": (("date", "row", "col"), [[[2.0, np.nan]]])}
        ),
    }

    lower, upper = validation._plot_limits(
        results,
        "dgimg",
        frame=0,
        symmetric=True,
    )

    assert lower == -upper
    assert lower < 0 < upper


def test_wic_variance_attachment_uses_half_count_offset():
    dataset = xr.Dataset(
        {
            "img": (
                ("date", "row", "col"),
                np.array([[[0.0, 10.0]], [[20.0, np.nan]]]),
            ),
            "id": xr.DataArray("WIC"),
        },
        coords={
            "date": [np.datetime64("2000-01-01"), np.datetime64("2000-01-02")],
            "row": [0],
            "col": [0, 1],
        },
    )

    result = validation._attach_measurement_variance(dataset, "WIC")
    coefficient = wic_variance_coefficient()
    expected = coefficient*dataset["img"] + 0.5*coefficient**2

    assert result["img_variance"].dims == dataset["img"].dims
    np.testing.assert_allclose(result["img_variance"], expected)
    assert result["img_variance"].attrs["poisson_offset"] == 0.5


def test_cli_parser_requires_explicit_sensor_paths_and_output():
    arguments = validation.build_parser().parse_args([
        "--wic", "wic.nc",
        "--si12", "si12.nc",
        "--si13", "si13.nc",
        "--output-dir", "output",
    ])

    assert arguments.wic == Path("wic.nc")
    assert arguments.si12 == Path("si12.nc")
    assert arguments.si13 == Path("si13.nc")
    assert arguments.output_dir == Path("output")
