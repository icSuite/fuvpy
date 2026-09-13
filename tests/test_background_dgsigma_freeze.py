import numpy as np
import xarray as xr

import fuvpy.src.background as background


def _synthetic_wic():
    """Return a small image sequence that requires multiple IRLS fits."""

    n_time = 3
    n_row = 6
    n_col = 8

    y, x = np.meshgrid(np.linspace(0, 1, n_row), np.linspace(0, 1, n_col), indexing="ij")
    sza_0 = 20 + 100*y
    dza_0 = 10 + 35*x
    sza = np.stack([sza_0 + frame for frame in range(n_time)])
    dza = np.stack([dza_0]*n_time)

    fraction = np.cos(np.deg2rad(sza))/np.cos(np.deg2rad(dza))
    time_change = 3*np.arange(n_time)[:, None, None]
    spatial_structure = 1.5*np.sin(5*x) + 0.8*np.cos(7*y)
    image = 150 + 45*fraction + 8*fraction**2 + time_change + spatial_structure

    dates = np.datetime64("2001-01-01") + 2*np.arange(n_time)*np.timedelta64(1, "m")
    return xr.Dataset(
        data_vars={
            "img": (("date", "row", "col"), image),
            "sza": (("date", "row", "col"), sza),
            "dza": (("date", "row", "col"), dza),
            "id": xr.DataArray("WIC")
        },
        coords={
            "date": dates,
            "row": np.arange(n_row),
            "col": np.arange(n_col)
        }
    )


def test_iteration_zero_dgsigma_is_retained_during_later_fits(monkeypatch):
    estimates = []

    def changing_noisemodel(fraction, image, model, weights, knots):
        estimate = 10.0 + len(estimates)
        estimates.append(estimate)
        return np.full_like(image, estimate)

    monkeypatch.setattr(background, "_noisemodel_binned", changing_noisemodel)

    result = background.backgroundmodel_BS(_synthetic_wic(), stop=1e-5)

    assert estimates == [10.0]
    np.testing.assert_allclose(result["dgsigma"].values[np.isfinite(result["dgsigma"])], 10.0)
    assert np.any(result["dgweight"].values < 1)


def test_dgsigma_freeze_iteration_is_not_configurable():
    with np.testing.assert_raises_regex(TypeError, "dgsigma_freeze_iteration"):
        background.backgroundmodel_BS(_synthetic_wic(), dgsigma_freeze_iteration=1)
