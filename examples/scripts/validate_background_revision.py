#!/usr/bin/env python3
"""Compare cumulative BS background revisions on complete native-grid orbits.

This script prepares figures and summary metrics for scientific review. Its
output is not evidence that the revised model has passed validation. Supply
complete orbit NetCDF files; short extracts do not reproduce the temporal fit.
"""

import argparse
import csv
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from fuvpy.src.background import backgroundmodel_BS
from fuvpy.src.measurement_uncertainty import (
    calibration_multiplier,
    si_measurement_variance,
)


DETECTOR_SHAPES = {
    "WIC": (256, 256),
    "SI12": (128, 128),
    "SI13": (128, 128),
}

VARIANTS = [
    ("released", {}),
    ("corrected_dgsigma", {"noise_model": "binned"}),
    (
        "bin_population",
        {
            "noise_model": "binned",
            "spatial_weighting": "bin_population",
        },
    ),
    (
        "frozen_dgsigma",
        {
            "noise_model": "binned",
            "spatial_weighting": "bin_population",
            "dgsigma_freeze_iteration": 0,
        },
    ),
    (
        "measurement_weighted",
        {
            "noise_model": "binned",
            "spatial_weighting": "bin_population",
            "dgsigma_freeze_iteration": 0,
            "input_variance": "img_variance",
            "variance_weighted": True,
        },
    ),
    (
        "monotonic",
        {
            "noise_model": "binned",
            "spatial_weighting": "bin_population",
            "dgsigma_freeze_iteration": 0,
            "input_variance": "img_variance",
            "variance_weighted": True,
            "monotonic": True,
            "monotonic_solver": "normal",
        },
    ),
]

PLOT_FIELDS = [
    ("dgmodel", "viridis", False),
    ("dgimg", "RdBu_r", True),
    ("dgsigma", "magma", False),
    ("dgweight", "viridis", False),
]


def _sensor_name(dataset):
    """Return one sensor identifier from a scalar or repeated id field."""

    values = np.asarray(dataset["id"]).ravel()
    names = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf8")
        names.append(str(value).strip().upper())

    unique = np.unique(names)
    if unique.size != 1:
        raise ValueError("Dataset must contain exactly one sensor id")

    return unique[0]


def _load_native_orbit(path, expected_sensor):
    """Load and check one complete native detector-grid orbit dataset."""

    dataset = xr.load_dataset(path)
    missing = [
        name for name in ["img", "date", "sza", "dza", "id"]
        if name not in dataset
    ]
    if missing:
        raise ValueError(f"{path} is missing required fields: {missing}")

    sensor = _sensor_name(dataset)
    if sensor != expected_sensor:
        raise ValueError(
            f"Expected {expected_sensor} dataset at {path}, found {sensor}"
        )

    required_dimensions = {"date", "row", "col"}
    for name in ["img", "sza", "dza"]:
        if set(dataset[name].dims) != required_dimensions:
            raise ValueError(
                f"{name} must have date, row, and col dimensions"
            )
        dataset[name] = dataset[name].transpose("date", "row", "col")

    detector_shape = (dataset.sizes["row"], dataset.sizes["col"])
    if detector_shape != DETECTOR_SHAPES[sensor]:
        raise ValueError(
            f"{sensor} must use its full native detector grid "
            f"{DETECTOR_SHAPES[sensor]}, found {detector_shape}"
        )
    if dataset.sizes["date"] < 2:
        raise ValueError("A complete orbit sequence must contain multiple frames")

    return dataset


def _attach_measurement_variance(dataset, sensor):
    """Preserve stored variance or attach the SI Step 4 variance model."""

    if "img_variance" in dataset:
        variance = dataset["img_variance"]
        if set(variance.dims) != {"date", "row", "col"}:
            raise ValueError(
                "Stored img_variance must have date, row, and col dimensions"
            )
        variance = variance.transpose("date", "row", "col")
        if variance.shape != dataset["img"].shape:
            raise ValueError("Stored img_variance must match the img shape")

        fraction = (
            np.cos(np.deg2rad(dataset["sza"]))
            / np.cos(np.deg2rad(dataset["dza"]))
        )
        fitted = (
            (dataset["sza"] >= 0)
            & np.isfinite(dataset["img"])
            & (fraction > -3.5)
            & (fraction < 3.5)
        )
        fitted_variance = np.asarray(variance)[np.asarray(fitted)]
        if np.any(~np.isfinite(fitted_variance)):
            raise ValueError("Stored img_variance must be finite at fitted pixels")
        if np.any(fitted_variance <= 0):
            raise ValueError(
                "Stored img_variance must be strictly positive at fitted pixels"
            )

        dataset["img_variance"] = variance
        return dataset

    if sensor == "WIC":
        raise ValueError(
            "WIC orbit lacks pre-reflattening img_variance; load and save the "
            "complete orbit with read_idl(..., measurement_variance=True, "
            "poisson_offset=0.5)"
        )

    image = np.asarray(dataset["img"])
    multiplier = calibration_multiplier(dataset["date"].values, sensor)
    variance = si_measurement_variance(
        image,
        multiplier,
        poisson_offset=0.5,
    )
    valid = np.isfinite(image) & np.isfinite(multiplier) & (multiplier > 0)
    dataset["calibration_multiplier"] = (
        ("date", "row", "col"),
        np.where(valid, multiplier, np.nan),
    )
    dataset["calibration_multiplier"].attrs = {
        "long_name": "FUVVIEW3 mission correction times flat field",
        "units": "1",
    }

    dataset["img_variance"] = (
        ("date", "row", "col"),
        np.where(valid, variance, np.nan),
    )
    dataset["img_variance"].attrs = {
        "long_name": "Conditional detector measurement variance of img",
        "units": "corrected count squared",
        "poisson_offset": 0.5,
        "excludes": (
            "calibration uncertainty, background-model discrepancy, and "
            "atmospheric-response uncertainty"
        ),
    }

    return dataset


def _representative_frames(n_frames, count):
    """Choose evenly spaced frame indices, including both endpoints."""

    count = min(count, n_frames)
    return np.unique(np.linspace(0, n_frames - 1, count, dtype=int))


def _plot_limits(results, field, frame, symmetric):
    """Return shared robust plot limits across all variants."""

    if field == "dgweight":
        return 0.0, 1.0

    values = np.concatenate([
        np.asarray(result[field].isel(date=frame)).ravel()
        for result in results.values()
    ])
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0

    lower, upper = np.percentile(values, [2, 98])
    if symmetric:
        limit = max(abs(lower), abs(upper))
        if limit == 0:
            limit = 1.0
        return -limit, limit
    if lower == upper:
        padding = max(abs(lower)*0.01, 1.0)
        return lower - padding, upper + padding

    return lower, upper


def _summary_metrics(sensor, variant, elapsed, result):
    """Return compact robust fit summaries for one variant."""

    corrected = np.asarray(result["dgimg"])
    corrected = corrected[np.isfinite(corrected)]
    corrected_median = np.median(corrected)

    weights = np.asarray(result["dgweight"])
    weights = weights[np.isfinite(weights)]

    return {
        "sensor": sensor,
        "variant": variant,
        "fit_seconds": elapsed,
        "dgmodel_median": float(np.nanmedian(result["dgmodel"])),
        "dgimg_median": float(corrected_median),
        "dgimg_mad": float(np.median(np.abs(corrected - corrected_median))),
        "dgsigma_median": float(np.nanmedian(result["dgsigma"])),
        "dgweight_median": float(np.median(weights)),
        "dgweight_below_half_fraction": float(np.mean(weights < 0.5)),
        "dgnorm_residual": float(result["dgnorm_residual"]),
        "dgnorm_model": float(result["dgnorm_model"]),
        "dgfit_iterations": int(result["dgfit_iterations"]),
        "dgfit_converged": bool(result["dgfit_converged"]),
    }


def _save_comparison_figures(sensor, results, output_directory, frame_count):
    """Save one compact four-field comparison for each selected frame."""

    first_result = next(iter(results.values()))
    frames = _representative_frames(first_result.sizes["date"], frame_count)

    for frame in frames:
        figure, axes = plt.subplots(
            len(VARIANTS),
            len(PLOT_FIELDS),
            figsize=(12, 14),
            constrained_layout=True,
        )

        for column, (field, colour_map, symmetric) in enumerate(PLOT_FIELDS):
            lower, upper = _plot_limits(
                results,
                field,
                frame,
                symmetric,
            )
            for row, (variant, _) in enumerate(VARIANTS):
                axis = axes[row, column]
                image = axis.imshow(
                    results[variant][field].isel(date=frame),
                    origin="lower",
                    cmap=colour_map,
                    vmin=lower,
                    vmax=upper,
                )
                axis.set_xticks([])
                axis.set_yticks([])
                if row == 0:
                    axis.set_title(field)
                if column == 0:
                    axis.set_ylabel(variant.replace("_", " "))
                if row == len(VARIANTS) - 1:
                    figure.colorbar(image, ax=axes[:, column], shrink=0.65)

        date = np.datetime_as_string(
            first_result["date"].values[frame],
            unit="s",
        )
        figure.suptitle(
            f"{sensor} frame {frame} ({date}) — comparison, not validation"
        )
        figure.savefig(
            output_directory / f"{sensor.lower()}_frame_{frame:04d}.png",
            dpi=150,
        )
        plt.close(figure)


def _write_metrics(rows, output_directory):
    """Write all sensor/variant summaries to one CSV file."""

    path = output_directory / "background_revision_metrics.csv"
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Compare cumulative BS revisions on complete native-grid WIC, "
            "SI12, and SI13 orbit NetCDF files. This prepares validation "
            "evidence; it does not declare the model validated."
        )
    )
    parser.add_argument("--wic", type=Path, required=True)
    parser.add_argument("--si12", type=Path, required=True)
    parser.add_argument("--si13", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=3)
    parser.add_argument("--damping", type=float, default=0.0)
    parser.add_argument("--tukey", type=float, default=5.0)
    parser.add_argument("--stop", type=float, default=1e-3)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.frames < 1:
        raise ValueError("--frames must be at least one")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_paths = {
        "WIC": args.wic,
        "SI12": args.si12,
        "SI13": args.si13,
    }
    common_settings = {
        "dampingVal": args.damping,
        "tukeyVal": args.tukey,
        "stop": args.stop,
        "diagnostics": True,
    }

    metric_rows = []
    for sensor, path in input_paths.items():
        dataset = _load_native_orbit(path, sensor)
        dataset = _attach_measurement_variance(dataset, sensor)

        results = {}
        for variant, settings in VARIANTS:
            started = perf_counter()
            result = backgroundmodel_BS(
                dataset,
                **common_settings,
                **settings,
            )
            elapsed = perf_counter() - started
            results[variant] = result
            metric_rows.append(
                _summary_metrics(sensor, variant, elapsed, result)
            )

        _save_comparison_figures(
            sensor,
            results,
            args.output_dir,
            args.frames,
        )

    _write_metrics(metric_rows, args.output_dir)
    print(
        "Comparison outputs written. Scientific validation has not been "
        "declared; review the complete-orbit evidence."
    )


if __name__ == "__main__":
    main()
