# Current State

Last reviewed: 2026-09-13
Branch: `background_model_asymmetry`
Status: BS-only implementation complete; production integration remains

## Active background architecture

fuvpy now has one background-removal function: `backgroundmodel_BS()`.
It supports:

- the ordinary one-dimensional BS background; and
- the selected directional model `B0(x) + T(SZA)*S*Bs(x)` through
  `spatial_model="legacy_x_azimuth"`.

The directional model uses `x=cos(SZA)/cos(DZA)`, the pole-safe solar-view
sine factor `S`, and the optional Spencer taper `T`. The cosine factor is
still calculated and saved for visualization. `monotonic=True` constrains the
base component, not the complete base-plus-directional model.

`backgroundmodel_SH()`, `backgroundmodel_SBS()`, their helpers and public
interfaces, and the rejected compressed-nightside `asymmetric_legacy`
experiment have been removed. New corrected-image workflows use `dgimg`.
Historical NetCDF files with SH fields remain readable as ordinary xarray
datasets.

## Candidate processing configuration

The current candidate combines:

```python
backgroundmodel_BS(
    images,
    spatial_model="legacy_x_azimuth",
    directional_x_knots=[-3.5, 0, 1.5, 3.5],
    spencer_taper=True,
    variance_weighted=True,
    monotonic=True
)
```

The production value of dimensionless `dampingVal` has not been selected.

## Retained revisions

- `dgsigma` is estimated from one statistic per populated geometry bin rather
  than repeating bin statistics per pixel. The obsolete legacy model has been
  removed.
- Every populated geometry bin receives equal total squared influence within
  a frame. The obsolete Gaussian spatial-weighting alternative has been
  removed.
- `dgsigma` is always estimated after the initial fit and then held fixed while
  later model and Tukey updates continue. This is no longer configurable.
- Tukey values enter as objective weights through `sqrt(Tukey)` row
  multipliers.
- BS damping is scaled by the median supported positive diagonal of the
  initial weighted normal matrix.
- `read_idl(..., measurement_variance=True)` attaches detector measurement
  variance before optional WIC reflattening.
- `variance_weighted=True` applies measurement precision to the BS inversion
  without including it in the empirical `dgsigma` calculation.
- The low-cost monotonic solver has explicit success/fallback flags, and every
  result stores scalar convergence diagnostics.

## Repository cleanup

The obsolete article examples, notebooks, paper figures, generated comparison
figures, Binder environment, empty submodule file, unused `image_orbit.h5`, and
generated caches were removed. Generated figures are now ignored.

The retained research support is:

- `scripts/compare_models.py`;
- `scripts/extract_legacy_background.py`;
- `scripts/validate_background_revision.py`;
- the WIC, SI12, and SI13 IDL validation sequences under `data/`; and
- compact historical BS and SH comparison fields under
  `data/legacy_background/`: `dgmodel`, `dgimg`, `dgweight`, `dgsigma`,
  `shmodel`, `shimg`, and `shweight`.

## Verification

The focused suite reports 38 passes. It was consolidated from 88 by combining
repeated error/parameter cases and removing tests of validation-script plotting
plumbing; package-level scientific behavior remains covered. `img_variance` now records what it was
calculated from, which uncertainty sources it excludes, and the deterministic
WIC reflattening assumption. No active source, test, or script refers to SH,
SBS, or `asymmetric_legacy`.

Monotonic solver failures no longer terminate processing. A successful
ordinary fallback is flagged separately from an all-NaN failed product, and
the saved monotonic metadata describes the result actually returned.

Every BS result now records its effective fit configuration and scalar
convergence state on `dgmodel`; solver outcome remains in `dgsolver_flag`.
Focused tests confirm that these values survive NetCDF save/reload.
`backgroundmodel_BS()` always fits `img`; the obsolete alternative `inImg`
path has been removed. It always returns a copy and rejects the obsolete
`inplace`, `diagnostics`, `input_variance`, `monotonic_solver`, `dzalim`, and
`minlat` arguments. Variance-weighted fits always read `img_variance`, and
production monotonic fits always use the low-cost normal-system solver.

## Production boundary

Do not start complete IMAGE processing yet:

1. Complete points 5--6 of
   `vault/03_Validation/FUVPY Production Reprocessing Checklist.md`.
2. Update icBuilder to call the revised loader and BS-only model and to use
   `dgimg`/`dgweight` instead of SH fields.
3. Make the orbit runner interruption-safe with atomic saving and output
   validation.

Detailed decisions, algorithms, and audit evidence remain in the dedicated
vault notes rather than this live status file.
