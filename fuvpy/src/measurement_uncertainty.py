"""Detector measurement variance for calibrated IMAGE/FUV images.

These calculations describe counting or detector-output noise conditional on
the supplied calibration. They do not include calibration uncertainty,
background-model discrepancy, or atmospheric-response uncertainty.
"""

from datetime import datetime, timedelta
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.io import readsav
from scipy.ndimage import rotate


# The published WIC equivalent-quanta value applies to one 2 x 2-pixel
# resolution cell, while the sensitivity is quoted per detector pixel.
WIC_SENSITIVITY_AD_PER_KR = 612.6
WIC_REFERENCE_BRIGHTNESS_KR = 0.1
WIC_PIXELS_PER_RESOLUTION_CELL = 4
WIC_EQUIVALENT_QUANTA = 23.0

CALIBRATION_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "data" / "calibration"
)

FLATFIELD_FILES = {
    "WIC": "wic_flatfield_dbase.idl",
    "SI12": "si12_flatfield_dbase.idl",
    "SI13": "si_flatfield_dbase.idl",
}

MISSION_FIELDS = {
    "WIC": "finalwic",
    "SI12": "finals12",
    "SI13": "finals13",
}

BOOM_LOSS_DATE = datetime(2000, 1, 1) + timedelta(days=277)


def _sensor_name(sensor):
    """Return a supported FUVVIEW3 sensor name."""

    sensor = str(sensor).strip().upper()
    if sensor not in FLATFIELD_FILES:
        raise ValueError(
            "Measurement variance is implemented for WIC, SI12, and SI13"
        )

    return sensor


def _poisson_offset(value):
    """Validate the optional count added to the plug-in mean."""

    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError("poisson_offset must be a non-negative number")

    if not np.isfinite(value) or value < 0:
        raise ValueError("poisson_offset must be a non-negative number")

    return value


@lru_cache(maxsize=None)
def _read_calibration_file(filename):
    """Read and cache one packaged FUVVIEW3 calibration file."""

    path = CALIBRATION_DIRECTORY / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing FUVVIEW3 calibration file: {path}")

    return readsav(str(path), python_dict=True)


def mission_correction(dates, sensor):
    """Return the date-dependent FUVVIEW3 mission-gain correction.

    The original routine silently returned one for dates absent from its
    table. Variance calculation instead raises because that would amount to
    inventing an unsupported calibration multiplier.
    """

    sensor = _sensor_name(sensor)
    dates = pd.DatetimeIndex(np.atleast_1d(dates))

    table = _read_calibration_file("fuv_mission_correction.sav")
    years = np.asarray(table["newyear"], dtype=int)
    days = np.asarray(table["newday"], dtype=int)
    sensitivity = np.asarray(table[MISSION_FIELDS[sensor]], dtype=float)

    reference = (years == 2000) & (days == 171)
    reference_sensitivity = sensitivity[reference][0]

    correction = np.empty(dates.size, dtype=float)
    for index, date in enumerate(dates):
        match = (years == date.year) & (days == date.dayofyear)
        if np.count_nonzero(match) != 1:
            raise ValueError(
                f"No unique {sensor} mission correction for {date.date()}"
            )

        correction[index] = reference_sensitivity / sensitivity[match][0]

    return correction


def _flatfield_dates(table):
    """Convert flat-field year/day columns to datetimes."""

    years = [int(value) for value in table["yearlist"]]
    days = [int(value) for value in table["daylist"]]

    return np.array([
        datetime(year, 1, 1) + timedelta(days=day - 1)
        for year, day in zip(years, days)
    ])


def flatfield_multiplier(date, sensor):
    """Reconstruct the FUVVIEW3 flat-field multiplier for one image.

    The nearest profile is selected without crossing the day-278 boom-loss
    boundary. Post-boom profiles receive the rotation and shift used for the
    final saved-image orientation.
    """

    sensor = _sensor_name(sensor)
    date = pd.Timestamp(date).normalize().to_pydatetime()

    table = _read_calibration_file(FLATFIELD_FILES[sensor])
    available_dates = _flatfield_dates(table)

    if date < BOOM_LOSS_DATE:
        available = available_dates < BOOM_LOSS_DATE
    else:
        available = available_dates >= BOOM_LOSS_DATE

    candidates = np.flatnonzero(available)
    distance = np.abs(available_dates[candidates] - date)
    selected = candidates[np.argmin(distance)]

    profile = np.asarray(table["flatfields"][:, selected], dtype=float)
    flatfield = np.repeat(profile[:, None], profile.size, axis=1)

    if date >= BOOM_LOSS_DATE:
        flatfield = rotate(
            flatfield,
            -1.85,
            reshape=False,
            order=1,
            mode="nearest",
            prefilter=False,
        )
        flatfield = np.roll(flatfield, 2, axis=0)

    return flatfield


def calibration_multiplier(dates, sensor):
    """Return mission correction times flat field for every image."""

    sensor = _sensor_name(sensor)
    dates = pd.DatetimeIndex(np.atleast_1d(dates))
    mission = mission_correction(dates, sensor)

    multiplier = []
    for date, mission_value in zip(dates, mission):
        multiplier.append(
            mission_value * flatfield_multiplier(date, sensor)
        )

    return np.asarray(multiplier)


def si_measurement_variance(image, multiplier, poisson_offset=0.0):
    """Return SI photon-counting variance in corrected-count squared.

    For saved image ``y = a*N`` and Poisson event count ``N``, the plug-in
    estimate is ``variance = a*y + poisson_offset*a**2``.
    """

    poisson_offset = _poisson_offset(poisson_offset)
    image = np.asarray(image)
    multiplier = np.asarray(multiplier)

    return multiplier*image + poisson_offset*multiplier**2


def wic_variance_coefficient():
    """Return the published reference WIC coefficient in AD units."""

    signal_per_pixel = (
        WIC_REFERENCE_BRIGHTNESS_KR * WIC_SENSITIVITY_AD_PER_KR
    )
    signal_per_resolution_cell = (
        WIC_PIXELS_PER_RESOLUTION_CELL * signal_per_pixel
    )

    return signal_per_resolution_cell / WIC_EQUIVALENT_QUANTA


def wic_measurement_variance(image, poisson_offset=0.0):
    """Return the reference WIC compound-Poisson variance in AD squared.

    The published equivalent-output-quanta relation gives
    ``variance = 10.65*image + poisson_offset*10.65**2`` (10.65 rounded).
    Additive CCD/read noise and sub-resolution-cell covariance are unknown and
    therefore excluded.
    """

    poisson_offset = _poisson_offset(poisson_offset)
    coefficient = wic_variance_coefficient()

    return coefficient*np.asarray(image) + poisson_offset*coefficient**2
