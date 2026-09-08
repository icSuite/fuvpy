# Handoff - Latest

Last updated: 2026-09-08
Branch: `background_model_revision`
Checkpoint: Steps 1 through 3 complete

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

## Verification

`pytest -q tests/test_background_dgsigma.py
tests/test_background_spatial_weighting.py
tests/test_background_dgsigma_freeze.py` passes 13 tests. The tests confirm
the Step 1 and 2 requirements, retain the iteration-0 scale while later robust
weights update, reproduce the released iteration-1 recalculation by default,
and reject invalid freeze settings.

## Next action

Implement Step 4 from `vault/02_Algorithm/Background Model Revision
Implementation Plan.md` in a separate commit. Do not add experimental
geometry, uncertainty, diagnostics, generated figures, or unrelated
abstractions.

## Portfolio impact

The clean background-model revision has completed its first three scientific
changes. This is an implementation checkpoint only; package defaults and
downstream production guidance are unchanged pending combined validation.
