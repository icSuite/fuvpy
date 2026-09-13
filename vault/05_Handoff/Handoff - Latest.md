# Handoff - Latest

Last updated: 2026-09-13
Branch: `background_model_asymmetry`
Checkpoint: BS-only architecture complete; numbered production checklist established

## What is implemented

The maintained background function is `backgroundmodel_BS()` with either the
ordinary one-dimensional model or the selected directional model:

```python
backgroundmodel_BS(
    images,
    spatial_model="legacy_x_azimuth",
    directional_x_knots=[-3.5, 0, 1.5, 3.5],
    spencer_taper=True,
    noise_model="binned",
    spatial_weighting="bin_population",
    dgsigma_freeze_iteration=0,
    input_variance="img_variance",
    variance_weighted=True,
    monotonic=True,
    diagnostics=True
)
```

The model is `B0(x) + T(SZA)*S*Bs(x)`. It saves the base and directional
components, both solar-view factors, and the Spencer taper. Monotonicity
applies only to the base component.

SH and SBS implementations, accessors, exports, helpers, and examples are
gone. The rejected compressed-nightside `asymmetric_legacy` path is also gone.
The corrected image produced by active processing is `dgimg`.

## Repository layout after cleanup

The article-specific `examples/` content and tracked generated figures were
removed because they describe an obsolete algorithm. The retained scripts are:

- `scripts/compare_models.py` — compare the historical BS fields, revised
  one-dimensional model, and selected directional model;
- `scripts/extract_legacy_background.py` — rebuild the compact historical
  comparison files from complete old orbit files; and
- `scripts/validate_background_revision.py` — compare cumulative BS revisions
  on complete native-grid orbits.

The IDL validation sequences and compact historical BS/SH comparison fields
remain under `data/`. Each compact file contains `dgmodel`, `dgimg`,
`dgweight`, `dgsigma`, `shmodel`, `shimg`, and `shweight`. Generated figures
are ignored. The obsolete Binder environment,
empty `.gitmodules`, unused 89 MB `image_orbit.h5`, and local build/test caches
were removed.

## Verification

`pytest -q` reports 93 passes. `img_variance` now records `calculated_from`,
`excludes`, and the WIC `reflattening_assumption`. The metadata fix does not
change the variance values or background fit. No active source, test, or
script refers to SH, SBS, or `asymmetric_legacy`.

`monotonic_fallback` now controls whether an ordinary fit is attempted after
a monotonic failure. Failed fits return flagged all-NaN products rather than
terminating processing. Flags distinguish successful fallback, no fallback,
and failed fallback, and the saved monotonic metadata reflects the returned
fit.

`dgmodel` now stores the effective BS configuration, including spline
breakpoints and orders, weighting and variance settings, monotonic settings,
damping value and scale, iteration count, and convergence state. The metadata
and solver flags survive NetCDF save/reload.

`backgroundmodel_BS()` now always fits `img`. The obsolete `inImg` option and
its redundant metadata were removed.

## Downstream incompatibility

icBuilder still calls `backgroundmodel_SH()` and expects `shimg` and
`shweight`. It cannot use this fuvpy branch until its orbit runner and detector
reader are changed to the BS-only fields.

## Next action

Use `vault/03_Validation/FUVPY Production Reprocessing Checklist.md` as the
numbered source of truth. Points 1--4 are complete. Point 5, the complete-orbit
pilot and damping selection, is next.

After points 5--6 are complete:

1. Update icBuilder to use the revised loader and BS-only background call.
2. Add atomic saving and restart validation to the orbit runner.

The production-readiness audit is in
`vault/03_Validation/Production Readiness Audit - 2026-09-12.md`.
