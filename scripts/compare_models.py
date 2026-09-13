"""Compare selected BS background-model options frame by frame.

Files separated by more than ten minutes are fitted as separate sequences so
that unrelated orbit segments do not share one temporal spline.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import fuvpy as fuv


REPOSITORY = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = REPOSITORY / "data"
LEGACY_BACKGROUND = EXAMPLE_DATA / "legacy_background"
OUTPUT_DIRECTORY = REPOSITORY / "figures" / "model_comparison"

SENSOR_DIRECTORIES = {
    "WIC": "wic_data",
    "SI12": "s12_data",
    "SI13": "s13_data",
}

BS_SETTINGS = {
    "stop": 0.01,
    "tukeyVal": 5,
    "dampingVal": 1e-2
}


def split_image_sequences(images, maximum_gap_minutes=10):
    """Split images when adjacent frames are from different orbit segments."""

    time_gap = np.diff(images.date.values) / np.timedelta64(1, "m")
    split_after = np.flatnonzero(time_gap > maximum_gap_minutes) + 1
    indices = np.split(np.arange(images.sizes["date"]), split_after)

    return [images.isel(date=index) for index in indices]


def colour_limits(values, symmetric=False):
    """Return readable limits without letting isolated pixels set the scale."""

    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0, 1

    if symmetric:
        limit = np.nanpercentile(np.abs(values), 99)
        return -limit, limit

    low, high = np.nanpercentile(values, [1, 99])
    if low == high:
        high = low + 1

    return low, high


def plot_image(axis, values, title, cmap, limits):
    """Plot one detector-grid field."""

    cmap = plt.get_cmap(cmap).copy()
    cmap.set_bad("black")

    axis.imshow(values, origin="lower", cmap=cmap, vmin=limits[0], vmax=limits[1])
    axis.set_title(title)
    axis.set_axis_off()


def plot_frame(original, legacy, m2, m3, frame, output_directory):
    """Plot all model variants for one frame."""
    
    # Img, Variance, SZA, DZA
    # x, cos, sine, spencer
    # org
    # new old
    # new
    # new : base + sine
    
    # x, cos, sine, spencer
    
    original_image = original.img.isel(date=frame).values
    variance = original.img_variance.isel(date=frame).values
    sza = original.sza.isel(date=frame).values
    dza = original.dza.isel(date=frame).values
    
    x = np.cos(np.deg2rad(original.sza.isel(date=frame))) / np.cos(np.deg2rad(original.dza.isel(date=frame)))
    cos = m3.solar_view_cosine_factor.isel(date=frame).values
    sin = m3.solar_view_sine_factor.isel(date=frame).values
    spencer = m3.spencer_taper_factor.isel(date=frame).values

    legacy_model = legacy.dgmodel.isel(date=frame).values
    legacy_model_sh = legacy.shmodel.isel(date=frame).values
    m2_model = m2.dgmodel.isel(date=frame).values
    m3_model = m3.dgmodel.isel(date=frame).values
    m3_model_base = m3.dgmodel_base.isel(date=frame).values
    m3_model_sin = m3.dgmodel_sine.isel(date=frame).values    
    model_values = np.stack([legacy_model, m2_model, m3_model, m3_model_base])
    model_asym_values = np.stack([legacy_model_sh, m3_model_sin])

    legacy_corrected = legacy.dgimg.isel(date=frame).values
    legacy_corrected_sh = legacy.shimg.isel(date=frame).values
    m2_corrected = m2.dgimg.isel(date=frame).values
    m3_corrected = m3.dgimg.isel(date=frame).values
    m3_corrected_base = m3.img.isel(date=frame).values - m3_model_base
    m3_corrected_sin = m3.img.isel(date=frame).values - m3_model_sin
    corrected_values = np.stack([legacy_corrected, m2_corrected, m3_corrected, m3_corrected_base, m3_corrected_sin, legacy_corrected_sh])

    legacy_weight = legacy.dgweight.isel(date=frame).values
    sh_weight = legacy.shweight.isel(date=frame).values
    m2_weight = m2.dgweight.isel(date=frame).values
    m3_weight = m3.dgweight.isel(date=frame).values

    legacy_uncertainty = legacy.dgsigma.isel(date=frame).values
    m2_uncertainty = m2.dgsigma.isel(date=frame).values
    m3_uncertainty = m3.dgsigma.isel(date=frame).values
    uncertainty_values = np.stack([legacy_uncertainty, m2_uncertainty, m3_uncertainty])

    image_limits = colour_limits(original_image)
    variance_limits = colour_limits(variance)
    model_limits = colour_limits(model_values, symmetric=True)
    asym_limits = colour_limits(model_asym_values, symmetric=True)
    corrected_limits = colour_limits(corrected_values, symmetric=True)
    uncertainty_limits = colour_limits(uncertainty_values, symmetric=True)

    figure, axes = plt.subplots(7, 4, figsize=(10, 17.5), sharex=True, sharey=True)

    panels = [
        (original_image, "Original image", "viridis", image_limits),
        (variance, "Estimated input variance", "magma", variance_limits),
        (sza, "SZA", "seismic", (0, 180)),
        (dza, "DZA", "magma", (0, 75)),
        (x, "cos(SZA)/cos(DZA)", "seismic", (-2, 2)),
        (cos, "Cosine factor", "seismic", (-1, 1)),
        (sin, "Sine factor", "seismic", (-1, 1)),
        (spencer, "Spencer taper", "magma", (0, 1)),
        (legacy_model, "Historical saved model", "seismic", model_limits),
        (legacy_corrected, "Corrected image", "seismic", corrected_limits),
        (legacy_weight, "Tukey weights", "viridis", (0, 1)),
        (legacy_uncertainty, "dgsigma", "seismic", uncertainty_limits),
        (legacy_model_sh, "Historical saved model (SH)", "seismic", asym_limits),
        (legacy_corrected_sh, "", "seismic", corrected_limits),
        (sh_weight, "", "viridis", (0, 1)),
        (legacy_uncertainty, "dgsigma (shared)", "seismic", uncertainty_limits),
        (m2_model, "Revised model", "seismic", model_limits),
        (m2_corrected, "", "seismic", corrected_limits),
        (m2_weight, "", "viridis", (0, 1)),
        (m2_uncertainty, "", "seismic", uncertainty_limits),
        (m3_model, "Asymmetric model", "seismic", model_limits),
        (m3_corrected, "", "seismic", corrected_limits),
        (m3_weight, "", "viridis", (0, 1)),
        (m3_uncertainty, "", "seismic", uncertainty_limits),
        (m3_model_base, "Base model", "seismic", model_limits),
        (m3_corrected_base, "Base corrected", "seismic", corrected_limits),
        (m3_model_sin, "Asymmetric model", "seismic", asym_limits),
        (m3_corrected_sin, "Asymmetric correction", "seismic", corrected_limits)        
    ]

    for axis, panel in zip(axes.flat, panels):
        plot_image(axis, *panel)

    sensor = str(original.id.values)
    date = original.date.isel(date=frame).values
    timestamp = np.datetime_as_string(date, unit="s")

    figure.suptitle(f"{sensor}  {timestamp}", y=.95)
    filename = timestamp.replace("-", "").replace(":", "") + ".png"
    figure.savefig(output_directory / filename, dpi=150, bbox_inches='tight')
    plt.close(figure)


def process_sensor(sensor, directory):
    """Load, fit, and plot every available frame for one sensor."""

    files = sorted((EXAMPLE_DATA / directory).glob("*.idl"))
    if not files:
        raise FileNotFoundError(f"No IDL files found for {sensor}")

    # Reflattening estimates one common WIC background from the loaded image
    # sequence. The example files contain partial orbits, so keep the original
    # FUVVIEW3 WIC flat field for this comparison.
    reflat = sensor == "WIC"
    #images = fuv.read_idl([str(path) for path in files], dzalim=75, reflat=reflat, measurement_variance=True, poisson_offset=0.5, wic_variance_model="correction_aware")
    images = fuv.read_idl([str(path) for path in files], dzalim=75, reflat=reflat, measurement_variance=True, poisson_offset=0.5, wic_variance_model='reference', viewing_geometry=True)

    legacy_file = LEGACY_BACKGROUND / f"{sensor.lower()}.nc"
    if not legacy_file.exists():
        raise FileNotFoundError(f"Missing {legacy_file}. Run extract_legacy_background.py first")
    with xr.open_dataset(legacy_file) as saved:
        legacy = saved.load()

    sequences = split_image_sequences(images)
    output_directory = OUTPUT_DIRECTORY / sensor.lower()
    output_directory.mkdir(parents=True, exist_ok=True)

    for sequence in sequences:
        try:
            legacy_sequence = legacy.sel(date=sequence.date)
        except KeyError as error:
            raise KeyError(f"Historical data do not contain every {sensor} frame in this sequence. Rerun extract_legacy_background.py") from error

        settings = dict(BS_SETTINGS)
        settings["n_tKnots"] = min(5, sequence.sizes["date"])
        settings["input_variance"] = "img_variance"
        settings["diagnostics"] = True
        settings['spatial_weighting'] = "bin_population"
        settings['dgsigma_freeze_iteration'] = 0

        # Combined revised model
        m2 = fuv.backgroundmodel_BS(sequence, noise_model="binned", variance_weighted=True, monotonic=True, **settings)
        
        # Asymmetric
        m3 = fuv.backgroundmodel_BS(sequence, spatial_model="legacy_x_azimuth", directional_x_knots=[-3.5, 0, 1.5, 3.5], spencer_taper=True, noise_model="binned", variance_weighted=True, monotonic=True, **settings)

        for frame in range(sequence.sizes["date"]):
            plot_frame(sequence, legacy_sequence, m2, m3, frame, output_directory)

    print(f"{sensor}: wrote {len(files)} figures to {output_directory}")


def main():
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    for sensor, directory in SENSOR_DIRECTORIES.items():
        process_sensor(sensor, directory)


if __name__ == "__main__":
    main()
