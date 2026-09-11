# Directional Background Model Implementation Plan

Last reviewed: 2026-09-11
Branch: `background_model_asymmetry`
Source of experimental evidence: `2D_model` at `1556ff1`
Status: Steps 1-6 implemented; representative complete-orbit comparison pending

## Goal

Add one restrained directional correction to the revised one-dimensional BS
background model without importing the experimental model zoo. The retained
candidate is

```text
B(x, t) = B0(x, t) + T(SZA) S Bs(x, t)
```

where

```text
x = cos(SZA)/cos(DZA)
S = solar_view_sine_factor
T = 1                           when spencer_taper=False
T = normalized Spencer profile when spencer_taper=True
```

`B0` remains the established legacy-coordinate model. `Bs` is a second
one-dimensional spline in the same coordinate. The signed geometric factor
`S` introduces a dawn-dusk distinction. There is no fitted cosine component.

The ordinary one-dimensional BS model must remain the default.

## Public interface

The intended model call is

```python
backgroundmodel_BS(images, spatial_model="legacy_x_azimuth", directional_x_knots=[-3.5, 0, 1.5, 3.5], spencer_taper=True)
```

The clean implementation exposes only three additions to the BS model:

- `spatial_model="legacy_x_azimuth"` selects the directional model;
- `directional_x_knots` controls the directional spline and defaults to the
  ordinary spatial knots;
- `spencer_taper=False/True` controls whether the Spencer profile multiplies
  the directional term.

Do not add an `include_cosine` switch. The cosine factor is a stored geometry
diagnostic, not a fitted term in this model.

## Implementation sequence

### 1. Recover and store the viewing geometry

Move the focused viewing-geometry calculation from `2D_model` into a small
module. `read_idl(..., viewing_geometry=True)` should use the FUVVIEW3 fields
`O_GCI`, `LAT_CTR`, `LON_CTR`, and `EMIS_HGT` together with pixel GLAT/GLON
and the subsolar point.

Calculate and save both of the pole-safe horizontal factors:

```text
C = horizontal_solar dot horizontal_viewing
S = local_up dot (horizontal_solar cross horizontal_viewing)
```

Store them as `solar_view_cosine_factor` and
`solar_view_sine_factor`. Both are required for visualization even though
only `S` enters the fitted model. Continue checking the reconstructed DZA
against the DZA stored by FUVVIEW3.

### 2. Add the Spencer taper

Copy the compact Spencer spherical-Chapman lookup table into package data and
add one interpolation helper. With `spencer_taper=False`, use `T=1`. With
`spencer_taper=True`, evaluate the normalized profile at each SZA and save it
as `spencer_taper_factor`.

The taper modifies only the directional term. It must not rescale `B0`.

### 3. Construct the minimal directional design

Keep the current legacy basis for `B0`. Construct a second one-dimensional
spline basis for `Bs` using `directional_x_knots`, multiply its rows by
`T*S`, and combine the two bases directly:

```text
G = [G0, (T*S) Gs]
```

Apply the existing temporal basis to both components. Keep this calculation
visible in `backgroundmodel_BS`; use a helper only for the distinct geometry
or spline-basis calculation.

### 4. Reuse the revised fitting machinery

The combined model should use the existing IRLS loop, binned `dgsigma`,
bin-population weighting, measurement-variance weighting, dimensionless
damping, convergence test, and fit-support mask. Do not introduce a second
fitting framework or the experimental base-first fitting sequence.

When `monotonic=True`, constrain only `B0`. The signed directional correction
means that the complete background is not globally monotonic in `x`; record
that limitation explicitly in the output metadata.

### 5. Save interpretable outputs

Retain the ordinary `dgmodel`, `dgimg`, `dgweight`, and `dgsigma` outputs. Add

- `dgmodel_base` for `B0`;
- `dgmodel_sine` for `T*S*Bs`;
- both solar-view factors from the loader; and
- `spencer_taper_factor` when the taper is enabled.

Record the selected spatial model, directional knots, taper state, and
monotonic constraint in attributes. Confirm that NetCDF writing and reading
preserve the new arrays and attributes.

### 6. Focused verification

Before processing complete orbits, add focused tests for:

1. reconstruction of stored DZA from bundled WIC, SI12, and SI13 files;
2. finite, correctly shaped cosine and sine factors;
3. the Spencer profile and both taper states;
4. exact recovery of a synthetic `B0 + T*S*Bs` background;
5. custom and default directional knots;
6. monotonic constraint acting on `B0` only;
7. compatibility with the revised weighting and uncertainty options;
8. unchanged output from the default one-dimensional model; and
9. NetCDF round-trip preservation.

### 7. Scientific comparison

After the implementation is stable, compare the ordinary and directional
models on representative complete orbits. Examine the fitted background,
corrected image, Tukey weights, base and directional components, geometry
factors, convergence, and sensitivity to the Spencer taper. Do not make the
directional model a default until that comparison is convincing.

## Explicitly excluded

Do not port the following experiments from `2D_model`:

- unsigned SZA-DZA tensor surfaces;
- signed-DZA and eastward-angle surfaces;
- fitted cosine components;
- SZA-only directional coefficients;
- two-dimensional sine-coordinate surfaces;
- full-replacement Spencer/sine or x/sine surfaces;
- scalar asymmetry amplitude;
- raised-cosine directional taper;
- base-only weight initialization and weight-recovery controls; or
- generic branching and validation machinery needed only by those models.

## Implementation checkpoint

Steps 1-6 were implemented on 2026-09-11. The clean implementation adds about
140 lines and removes 17 lines in `background.py` relative to
`background_model_revision`, compared with more than 800 added lines in the
experimental model-zoo version.

Twenty-one focused directional, geometry, taper, compatibility, and NetCDF
tests pass. The complete repository suite reports 100 passes and the same
three failures already recorded on `background_model_revision`: one deferred
asymmetric-coordinate knot test and two missing measurement-variance metadata
attributes. No complete IMAGE orbit was processed.
