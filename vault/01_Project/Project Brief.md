# Project Brief

Last reviewed: 2026-09-01

## Purpose

`fuvpy` is a Python package for global far-ultraviolet auroral images. It reads
IMAGE/FUV and Polar/UVI products into xarray datasets, fits and removes
background emissions, models auroral boundaries, and provides plotting and
inspection tools.

The `background_uncertainty` branch has one additional purpose: estimate the
uncertainty introduced by the fitted background removal so downstream IMAGE
products can propagate it through sensor coregistration, precipitation, and
conductance.

## Scope

The repository contains:

- file loading and WIC reflattening in `fuvpy/src/load.py`;
- dayglow and residual-background models in `fuvpy/src/background.py`;
- auroral-boundary models in `fuvpy/src/boundaries.py`;
- xarray accessors and plotting helpers;
- utilities for spherical harmonics, sunlight, coordinates, dates, and oval
  models; and
- notebooks, scripts, figures, and a small WIC example sequence.

The uncertainty work is limited initially to the existing BS and SH background
fits. It does not redesign the mean background model, calibrate the cameras,
or choose the production electron-energy method in icBuilder.

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
declared in `pyproject.toml`. The supplied environment is pinned to Python
3.10-era scientific packages. No automated tests are present at the initial
checkpoint.

## External ownership

The local repository tracks the public `aohma/fuvpy` upstream. The uncertainty
branch is local and has not been pushed. Do not assume permission to publish
changes upstream or change the remote.
