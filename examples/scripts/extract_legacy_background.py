"""Extract the historical background fields for the bundled IDL examples.

The complete orbit files are too large for the repository. This script keeps
only the dates represented by the example IDL files and the four fields used
by compare_models.py.
"""

from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr


REPOSITORY = Path(__file__).resolve().parents[2]
SOURCE_DIRECTORY = Path.home() / "Downloads" / "data"
OUTPUT_DIRECTORY = REPOSITORY / "data" / "legacy_background"

SENSORS = {
    "wic": ("wic", "wic_data"),
    "si12": ("s12", "s12_data"),
    "si13": ("s13", "s13_data")
}

FIELDS = ["dgmodel", "dgimg", "dgweight", "dgsigma"]


def idl_minute(path):
    """Read YYYY-DDD-HH-MM from the end of an IMAGE IDL filename."""

    timestamp = datetime.strptime(path.stem[-11:], "%Y%j%H%M")
    return np.datetime64(timestamp, "m")


def extract_sensor(sensor, source_prefix, idl_directory):
    """Extract the historical fields for one sensor."""

    idl_files = sorted((REPOSITORY / "data" / idl_directory).glob("*.idl"))
    wanted_minutes = np.array([idl_minute(path) for path in idl_files])
    extracted = []

    for source_file in sorted(SOURCE_DIRECTORY.glob(f"{source_prefix}_or*.nc")):
        with xr.open_dataset(source_file) as source:
            source_minutes = source.date.values.astype("datetime64[m]")
            keep = np.flatnonzero(np.isin(source_minutes, wanted_minutes))
            if keep.size == 0:
                continue

            orbit = int(source_file.stem[-4:])
            part = source[FIELDS].isel(date=keep).astype("float32").load()
            part["orbit"] = ("date", np.full(keep.size, orbit, dtype="int16"))
            extracted.append(part)

    if not extracted:
        print(f"{sensor}: no complete source files found")
        return

    result = xr.concat(extracted, dim="date").sortby("date")
    found_minutes = result.date.values.astype("datetime64[m]")
    missing_minutes = np.setdiff1d(wanted_minutes, found_minutes)

    result.attrs["description"] = "Historical fuvpy background results for the bundled IDL examples"
    result.attrs["source"] = str(SOURCE_DIRECTORY)
    result.attrs["selection"] = "Frames matched to example IDL filenames at one-minute precision"

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIRECTORY / f"{sensor}.nc"
    partial_file = output_file.with_suffix(".partial.nc")
    encoding = {field: {"zlib": True, "complevel": 6, "shuffle": True} for field in FIELDS}
    result.to_netcdf(partial_file, encoding=encoding)
    partial_file.replace(output_file)

    print(f"{sensor}: saved {result.sizes['date']} of {len(wanted_minutes)} frames to {output_file}")
    if missing_minutes.size:
        print(f"{sensor}: {missing_minutes.size} frames are still missing; rerun after the remaining orbit files arrive")


def main():
    for sensor, (source_prefix, idl_directory) in SENSORS.items():
        extract_sensor(sensor, source_prefix, idl_directory)


if __name__ == "__main__":
    main()
