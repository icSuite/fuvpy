# Current State

Last reviewed: 2026-09-08
Branch: `background_model_revision`
Implementation status: Steps 1 through 5 complete

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
- `dgsigma_freeze_iteration` selects the last zero-based IRLS iteration that
  may replace the operational residual-spread field. Setting it to zero keeps
  the initial-fit estimate fixed while later model and Tukey-weight updates
  continue normally. The default remains one to reproduce released behavior.
- `read_idl(..., measurement_variance=True)` adds conditional detector
  measurement variance for WIC, SI12, and SI13. Ordinary loading is unchanged.
- SI uses `a*img + poisson_offset*a**2`, where `a` is reconstructed from the
  packaged FUVVIEW3 mission-gain and flat-field tables. Reference WIC uses
  `10.65*img + poisson_offset*10.65**2` (coefficient rounded in prose).
- Variance is calculated before optional WIC reflattening and retained under
  the conditional assumption that its estimated additive background is fixed.
  Calibration, atmospheric-response, and background-model uncertainties are
  excluded explicitly in the metadata.
- `backgroundmodel_BS(..., input_variance="img_variance",
  variance_weighted=True)` applies `1/sqrt(variance)` only to the background
  inversion. Robust Tukey, spatial, and measurement row multipliers remain
  separate and visible. Empirical `dgsigma` and Tukey residual scaling do not
  use measurement variance.
- Measurement variance must be finite and strictly positive at every fitted
  pixel. Nonzero damping must be reassessed because absolute measurement
  weights alter the data-to-regularization balance. The released unweighted
  fit remains the default.

## Verification

The focused tests through Step 5 pass 35 tests. They cover the
dgsigma bin structure, known variance relation, and population invariance,
plus equal within-bin weights, equal bin influence, frame independence, and
excluded/empty-bin behavior, the legacy-coordinate scope, iteration-0 scale
retention, the released iteration-1 default, detector variance formulas,
calibration-table selection, loader opt-in behavior, masks and metadata, and
unsupported calibration states, plus constant-variance equivalence, reduced
high-variance influence, dgsigma separation, and invalid fit variance.

## Current boundary and next action

No geometry, background covariance, diagnostics, or monotonic-solver changes
are present on this branch. The next planned change is Step 6, adding the
optional monotonic spatial solver, as a separate commit.
