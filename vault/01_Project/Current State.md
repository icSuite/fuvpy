# Current State

Last reviewed: 2026-09-08
Branch: `background_model_revision`
Implementation status: Step 1 complete

## Confirmed implementation

- The released BS background calculation remains the default.
- `backgroundmodel_BS(..., noise_model="binned")` enables the corrected
  empirical `dgsigma` regression while retaining the legacy coordinate
  `cos(SZA)/cos(DZA)`.
- The corrected regression uses one nightside bin below zero and 14 dayside
  bins from 0 to 3.5 at width 0.25.
- Each populated bin contributes one weighted mean model-squared value and one
  weighted mean residual-squared value. These bin-level pairs are fitted as
  `residual_squared = a*model_squared + b` with nonnegative `a` and `b`.
- The fitted spread is evaluated at every fitted pixel as
  `sqrt(a*dgmodel**2 + b)`. Bin statistics are not expanded per pixel before
  regression.

## Verification

The focused `tests/test_background_dgsigma.py` checks the bin structure,
recovery of a known nonnegative variance relation, and invariance when the
population of one bin is duplicated. All three tests pass.

## Current boundary and next action

No geometry, spatial weighting, freeze-iteration, measurement-variance,
uncertainty, diagnostics, or monotonic-solver changes are present on this
branch. The next planned change is Step 2, bin-population spatial weighting,
as a separate commit.
