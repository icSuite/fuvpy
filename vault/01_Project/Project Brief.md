# Project Brief

Last reviewed: 2026-09-12

## Purpose

`fuvpy` is a Python package for global far-ultraviolet auroral images. It reads
IMAGE/FUV and Polar/UVI products into xarray datasets, fits and removes
background emissions, models auroral boundaries, and provides plotting and
inspection tools.

The `background_model_asymmetry` branch develops the revised BS background
model and its optional directional correction for downstream IMAGE products.

## Scope

The repository contains:

- file loading and WIC reflattening in `fuvpy/src/load.py`;
- the B-spline background model in `fuvpy/src/background.py`;
- auroral-boundary models in `fuvpy/src/boundaries.py`;
- xarray accessors and plotting helpers;
- utilities for sunlight, coordinates, dates, and oval models; and
- validation scripts, compact historical results, and three IDL validation
  sequences.

The background work is limited to `backgroundmodel_BS`, with or without its
directional asymmetry. It does not calibrate the cameras or choose the
production electron-energy method in icBuilder.

## Downstream role

icBuilder uses fuvpy-derived WIC, SI12, and SI13 detector images. A defensible
background-fit uncertainty could become one component of detector-level
Product-1 uncertainty and could be propagated by Monte Carlo into the
conductance ensemble used by icAnalyzer.

This component does not include detector counting noise, calibration and
flat-field uncertainty, atmospheric-response uncertainty, or background-model
discrepancy unless those sources are added explicitly.

## Dependencies and environment

`pyproject.toml` declares matplotlib, netCDF4, NumPy, pandas, SciPy, and xarray.
The README also requires `polplot`, which is imported by the source but is not
declared in `pyproject.toml`. Focused automated tests cover the revised
background model, loading, measurement variance, and viewing geometry.

## External ownership

`origin` is `icSuite/fuvpy`; the original `aohma/fuvpy` repository is retained
as `upstream`. The active branch is `background_model_asymmetry`.
