import numpy as np
import pytest
import xarray as xr

import fuvpy.src.load as load
from fuvpy.src.measurement_uncertainty import (
    calibration_multiplier,
    mission_correction,
    si_measurement_variance,
    wic_measurement_variance,
    wic_variance_coefficient,
)


def test_si_measurement_variance_and_low_count_offset():
    image = np.array([10.0, 0.0])
    multiplier = np.array([1.18, 1.18])

    variance = si_measurement_variance(
        image,
        multiplier,
        poisson_offset=0.5,
    )

    np.testing.assert_allclose(
        variance,
        multiplier*image + 0.5*multiplier**2,
    )


def test_wic_reference_variance_and_low_count_offset():
    image = np.array([10.0, 0.0])
    coefficient = 4*0.1*612.6/23

    variance = wic_measurement_variance(image, poisson_offset=0.5)

    np.testing.assert_allclose(wic_variance_coefficient(), coefficient)
    np.testing.assert_allclose(
        variance,
        coefficient*image + 0.5*coefficient**2,
    )


def test_packaged_tables_select_representative_mission_multipliers():
    for sensor in ["WIC", "SI12", "SI13"]:
        np.testing.assert_allclose(
            mission_correction(["2000-06-19"], sensor),
            [1.0],
        )

    np.testing.assert_allclose(
        mission_correction(["2000-05-16"], "SI12"),
        [0.9908789386401327],
    )
    multiplier = calibration_multiplier(["2000-05-16"], "SI12")
    assert multiplier.shape == (1, 128, 128)
    np.testing.assert_allclose(multiplier[0, :, 0], multiplier[0, :, -1])


@pytest.mark.parametrize("value", [-1, np.inf, "invalid"])
def test_poisson_offset_must_be_nonnegative_and_finite(value):
    with pytest.raises(ValueError, match="non-negative number"):
        si_measurement_variance([1.0], [1.0], poisson_offset=value)


def _measurement_dataset(sensor="SI12"):
    return xr.Dataset(
        data_vars={
            "img": (("date", "row", "col"), [[[10.0, 0.0], [20.0, 30.0]]]),
            "mlat": (("date", "row", "col"), [[[60.0, 61.0], [np.nan, 63.0]]]),
            "id": xr.DataArray(sensor),
        },
        coords={
            "date": [np.datetime64("2000-05-16")],
            "row": [0, 1],
            "col": [0, 1],
        },
    )


def test_si_variance_follows_image_dimensions_mask_and_metadata(monkeypatch):
    multiplier = np.array([[[1.0, 2.0], [3.0, 4.0]]])
    monkeypatch.setattr(
        load,
        "calibration_multiplier",
        lambda dates, sensor: multiplier,
    )

    result = load._add_measurement_variance(
        _measurement_dataset(),
        calibration_flags=[1],
        airglow_scales=[[0, 0, 0]],
        poisson_offset=0.5,
    )

    expected = multiplier*result["img"].values + 0.5*multiplier**2
    expected[0, 1, 0] = np.nan

    assert result["img_variance"].dims == result["img"].dims
    np.testing.assert_allclose(result["img_variance"], expected)
    assert np.isnan(result["calibration_multiplier"][0, 1, 0])
    assert result["img_variance"].attrs["poisson_offset"] == 0.5
    assert "calibration uncertainty" in result["img_variance"].attrs["excludes"]


@pytest.mark.parametrize(
    "sensor, flag, airglow, message",
    [
        ("SI12", 3, [0, 0, 0], "calibration flag 1"),
        ("WIC", 2, [0, 0, 0], "calibration flag 1 or 3"),
        ("SI13", 1, [1, 0, 0], "without FUVVIEW3 airglow"),
        ("VIS", 1, [0, 0, 0], "implemented for WIC, SI12, and SI13"),
    ],
)
def test_unsupported_calibration_states_are_rejected(
    sensor,
    flag,
    airglow,
    message,
):
    with pytest.raises(ValueError, match=message):
        load._add_measurement_variance(
            _measurement_dataset(sensor),
            calibration_flags=[flag],
            airglow_scales=[airglow],
            poisson_offset=0,
        )


def _mock_imageinfo(include_calibration=True):
    shape = (2, 2)
    imageinfo = {
        "inst_id": np.array([b"WIC"]),
        "image": np.array([np.full(shape, 10.0)]),
        "mlat": np.array([np.full(shape, 70.0)]),
        "mlon": np.array([np.zeros(shape)]),
        "mlt": np.array([np.zeros(shape)]),
        "glat": np.array([np.full(shape, 60.0)]),
        "glon": np.array([np.zeros(shape)]),
        "dza": np.array([np.full(shape, 30.0)]),
        "sza": np.array([np.full(shape, 100.0)]),
        "time": np.array([[2000137, 0]]),
    }
    if include_calibration:
        imageinfo["calibration_flag"] = np.array([3])
        imageinfo["airglow_scale"] = np.array([[0.0, 0.0, 0.0]])

    return imageinfo


def test_read_idl_variance_precedes_optional_wic_reflattening(monkeypatch):
    monkeypatch.setattr(
        load,
        "readsav",
        lambda filename: {"imageinfo": _mock_imageinfo()},
    )
    monkeypatch.setattr(
        load,
        "calibration_multiplier",
        lambda dates, sensor: np.ones((1, 2, 2)),
    )

    def add_known_reflattening(imgs):
        imgs["img"] = imgs["img"] + 100
        return imgs

    monkeypatch.setattr(load, "_reflatWIC", add_known_reflattening)

    result = load.read_idl(
        "mock.idl",
        remove_bad=False,
        measurement_variance=True,
        poisson_offset=0.5,
    )

    coefficient = wic_variance_coefficient()
    np.testing.assert_allclose(result["img"], 110.0)
    np.testing.assert_allclose(
        result["img_variance"],
        coefficient*10.0 + 0.5*coefficient**2,
    )
    assert int(result["calibration_flag"][0]) == 3
    assert "deterministic" in result["img_variance"].attrs[
        "reflattening_assumption"
    ]


def test_read_idl_default_does_not_require_calibration_fields(monkeypatch):
    monkeypatch.setattr(
        load,
        "readsav",
        lambda filename: {
            "imageinfo": _mock_imageinfo(include_calibration=False)
        },
    )

    result = load.read_idl("mock.idl", remove_bad=False, reflat=False)

    assert "img_variance" not in result
    assert "calibration_multiplier" not in result
