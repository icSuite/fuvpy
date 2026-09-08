# Handoff - Latest

Last updated: 2026-09-08
Branch: `background_model_revision`
Checkpoint: Steps 1 through 4 complete

## Completed

The released `backgroundmodel_BS` now has an explicit
`noise_model="binned"` option. It retains the legacy dayglow coordinate and
all other released fit behavior, but corrects the empirical `dgsigma`
regression by fitting one point per populated geometry bin. The default
`noise_model="legacy"` preserves the released result for comparison during
combined validation.

The bin edges define one fitted nightside bin below zero and 14 dayside bins
from 0 to 3.5 in steps of 0.25. Within each populated bin the calculation
forms the weighted means of `dgmodel**2` and `(img-dgmodel)**2`. A bounded
linear fit constrains both slope and intercept to be nonnegative, and
`dgsigma` is the square root of the fitted variance at every pixel.

The BS fit also has an explicit `spatial_weighting="bin_population"` option.
It bins the legacy coordinate at width 0.1 independently in each frame and
assigns each fitted row a multiplier proportional to `1/sqrt(N_bin)`. The
per-frame multipliers are normalized to mean square one over fitted pixels.
The released inverse-Gaussian multiplier remains the default through
`spatial_weighting="gaussian"`. The new option is explicitly rejected for VIS,
which does not use the legacy FUV coordinate. The selected spatial multiplier
continues to enter the `dgsigma` calculation through the same combined
robust/spatial weight used by the released path.

`dgsigma_freeze_iteration` now controls the last zero-based IRLS iteration
that recalculates the operational scale. With zero, the scale is estimated
after the initial fit and retained while subsequent model and Tukey-weight
updates proceed. The default is one, preserving the released
`iteration <= 1` sequence. The setting must be a nonnegative integer; no
diagnostic history machinery was added.

`read_idl(..., measurement_variance=True)` now reconstructs the original
FUVVIEW3 calibration multiplier and adds `img_variance` for WIC, SI12, and
SI13. SI uses `a*img + poisson_offset*a**2`; WIC uses the published reference
compound-Poisson coefficient derived from `4*0.1*612.6/23`, approximately
10.65 AD. The optional offset defaults to zero and is recorded in attributes.

The loader accepts WIC calibration flags 1 and 3 under the documented
assumption that `IMAGE` is the flag-1-equivalent corrected-AD field and flag 3
labels the separate `CIMAGE` product. SI requires flag 1, all sensors require
zero FUVVIEW3 airglow subtraction, and unsupported sensors or states raise.
Variance is calculated before optional WIC reflattening and kept unchanged
under the conditional deterministic-additive interpretation. Package data now
include the original small mission and flat-field tables.

## Verification

`pytest -q tests/test_background_dgsigma.py
tests/test_background_spatial_weighting.py
tests/test_background_dgsigma_freeze.py
tests/test_measurement_uncertainty.py` passes 26 tests. In addition to Steps
1--3, tests cover SI and WIC equations and offsets, representative mission and
flat-field selection, dimensions, masks, metadata, rejected calibration
states, pre-reflattening calculation, and unchanged default loader behavior.

## Next action

Implement Step 5 from `vault/02_Algorithm/Background Model Revision
Implementation Plan.md` in a separate commit. Do not add experimental
geometry, uncertainty, diagnostics, generated figures, or unrelated
abstractions.

## Portfolio impact

The clean background-model revision has completed its first four scientific
changes. This is an implementation checkpoint only; package defaults and
downstream production guidance are unchanged pending combined validation.
