# Current State

Last reviewed: 2026-09-08
Branch: `background_model_revision`
Implementation status: Steps 1 and 2 complete

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
- `backgroundmodel_BS(..., spatial_weighting="bin_population")` enables
  per-frame spatial sampling weights while the released Gaussian weighting
  remains the default. The option is limited to sensors using the legacy FUV
  coordinate; VIS is rejected explicitly.
- The legacy coordinate is divided into 0.1-wide bins. A fitted pixel in a
  bin of population `N` receives a row multiplier proportional to
  `1/sqrt(N)`. Per-frame normalization makes the mean squared multiplier one,
  so every populated bin has equal total squared influence without changing
  the overall data-to-damping scale.

## Verification

The focused dgsigma and spatial-weighting tests pass 7 tests. They cover the
dgsigma bin structure, known variance relation, and population invariance,
plus equal within-bin weights, equal bin influence, frame independence, and
excluded/empty-bin behavior, and the legacy-coordinate scope.

## Current boundary and next action

No geometry, freeze-iteration, measurement-variance, uncertainty, diagnostics,
or monotonic-solver changes are present on this branch. The next planned
change is Step 3, freezing `dgsigma` after iteration 0, as a separate commit.
