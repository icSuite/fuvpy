"""Compare selected BS background-model options frame by frame.

Files separated by more than ten minutes are fitted as separate sequences so
that unrelated orbit segments do not share one temporal spline.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import fuvpy as fuv


REPOSITORY = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = REPOSITORY / "data"
OUTPUT_DIRECTORY = REPOSITORY / "figures" / "bs_model_comparison"

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


def plot_frame(original, m1, m2, m3, m4, m5, m6, frame, output_directory):
    """Plot all model variants for one frame."""

    original_image = original.img.isel(date=frame).values
    sza = original.sza.isel(date=frame).values
    coordinate = np.cos(np.deg2rad(original.sza.isel(date=frame))) / np.cos(np.deg2rad(original.dza.isel(date=frame)))
    variance = original.img_variance.isel(date=frame).values

    m1_model = m1.dgmodel.isel(date=frame).values
    m2_model = m2.dgmodel.isel(date=frame).values
    m3_model = m3.dgmodel.isel(date=frame).values
    m4_model = m4.dgmodel.isel(date=frame).values
    m5_model = m5.dgmodel.isel(date=frame).values
    m6_model = m6.dgmodel.isel(date=frame).values
    model_values = np.stack([m1_model, m2_model, m3_model, m4_model, m5_model, m6_model])

    m1_corrected = m1.dgimg.isel(date=frame).values
    m2_corrected = m2.dgimg.isel(date=frame).values
    m3_corrected = m3.dgimg.isel(date=frame).values
    m4_corrected = m4.dgimg.isel(date=frame).values
    m5_corrected = m5.dgimg.isel(date=frame).values
    m6_corrected = m6.dgimg.isel(date=frame).values
    corrected_values = np.stack([m1_corrected, m2_corrected, m3_corrected, m4_corrected, m5_corrected, m6_corrected])

    m1_weight = m1.dgweight.isel(date=frame).values
    m2_weight = m2.dgweight.isel(date=frame).values
    m3_weight = m3.dgweight.isel(date=frame).values
    m4_weight = m4.dgweight.isel(date=frame).values
    m5_weight = m5.dgweight.isel(date=frame).values
    m6_weight = m6.dgweight.isel(date=frame).values

    m1_uncertainty = m1.dgsigma.isel(date=frame).values
    m2_uncertainty = m2.dgsigma.isel(date=frame).values
    m3_uncertainty = m3.dgsigma.isel(date=frame).values
    m4_uncertainty = m4.dgsigma.isel(date=frame).values
    m5_uncertainty = m5.dgsigma.isel(date=frame).values
    m6_uncertainty = m6.dgsigma.isel(date=frame).values
    uncertainty_values = np.stack([m1_uncertainty, m2_uncertainty, m3_uncertainty, m4_uncertainty, m5_uncertainty, m6_uncertainty])

    image_limits = colour_limits(original_image)
    variance_limits = colour_limits(variance)
    model_limits = colour_limits(model_values)
    corrected_limits = colour_limits(corrected_values, symmetric=True)
    uncertainty_limits = colour_limits(uncertainty_values, symmetric=True)

    figure, axes = plt.subplots(7, 4, figsize=(10, 17), sharex=True, sharey=True)

    panels = [
        (original_image, "Original corrected image", "viridis", image_limits),
        (variance, "Estimated input variance", "magma", variance_limits),
        (sza, "SZA", "magma", (0, 150)),
        (coordinate, "cos(SZA)/cos(DZA)", "seismic", (-2, 2)),
        (m1_model, "BS model : Original", "viridis", model_limits),
        (m1_corrected, "Corrected image", "seismic", corrected_limits),
        (m1_weight, "Tukey weights", "viridis", (0, 1)),
        (m1_uncertainty, "dgsigma", "seismic", uncertainty_limits),
        (m2_model, "Original + dgsigma", "viridis", model_limits),
        (m2_corrected, "", "seismic", corrected_limits),
        (m2_weight, "", "viridis", (0, 1)),
        (m2_uncertainty, "", "seismic", uncertainty_limits),
        (m3_model, "Original + variance", "viridis", model_limits),
        (m3_corrected, "", "seismic", corrected_limits),
        (m3_weight, "", "viridis", (0, 1)),
        (m3_uncertainty, "", "seismic", uncertainty_limits),
        (m4_model, "Original + mono", "viridis", model_limits),
        (m4_corrected, "", "seismic", corrected_limits),
        (m4_weight, "", "viridis", (0, 1)),
        (m4_uncertainty, "", "seismic", uncertainty_limits),
        (m5_model, "New", "viridis", model_limits),
        (m5_corrected, "", "seismic", corrected_limits),
        (m5_weight, "", "viridis", (0, 1)),
        (m5_uncertainty, "", "seismic", uncertainty_limits),
        (m6_model, "New : asymmetric", "viridis", model_limits),
        (m6_corrected, "", "seismic", corrected_limits),
        (m6_weight, "", "viridis", (0, 1)),
        (m6_uncertainty, "", "seismic", uncertainty_limits),
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
    images = fuv.read_idl([str(path) for path in files], dzalim=75, reflat=reflat, measurement_variance=True, poisson_offset=0.5, wic_variance_model='reference')

    sequences = split_image_sequences(images)
    output_directory = OUTPUT_DIRECTORY / sensor.lower()
    output_directory.mkdir(parents=True, exist_ok=True)

    for sequence in sequences:
        settings = dict(BS_SETTINGS)
        settings["n_tKnots"] = min(5, sequence.sizes["date"])
        settings["input_variance"] = "img_variance"
        settings["diagnostics"] = True
        settings['spatial_weighting'] = "bin_population"
        settings['dgsigma_freeze_iteration'] = 0

        # Original
        m1 = fuv.backgroundmodel_BS(sequence, **settings)

        # Original + corrected dgsigma
        m2 = fuv.backgroundmodel_BS(sequence, noise_model="binned", **settings)

        # Original + variance
        m3 = fuv.backgroundmodel_BS(sequence, noise_model="binned", variance_weighted=True, **settings)

        # Original + monotonic
        m4 = fuv.backgroundmodel_BS(sequence, noise_model="binned", monotonic=True, **settings)

        # Combined revised model
        m5 = fuv.backgroundmodel_BS(sequence, noise_model="binned", variance_weighted=True, monotonic=True, **settings)

        # Combined revised model with compressed nightside coordinates
        m6 = fuv.backgroundmodel_BS(sequence, noise_model="binned", variance_weighted=True, monotonic=True, asymmetric_legacy=True, **settings)

        for frame in range(sequence.sizes["date"]):
            plot_frame(sequence, m1, m2, m3, m4, m5, m6, frame, output_directory)

    print(f"{sensor}: wrote {len(files)} figures to {output_directory}")


def main():
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    for sensor, directory in SENSOR_DIRECTORIES.items():
        process_sensor(sensor, directory)


if __name__ == "__main__":
    main()
