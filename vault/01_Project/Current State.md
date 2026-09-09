# Current State

Last reviewed: 2026-09-09
Branch: `background_model_revision`
Implementation status: Steps 1 through 6 complete
Validation status: combined complete-orbit script ready; execution blocked

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
- `backgroundmodel_BS(..., asymmetric_legacy=True)` optionally preserves
  nonnegative `cos(SZA)/cos(DZA)` values and maps negative values to
  `nightside_width*tanh(x/nightside_width)`. The default width is 0.2; it must
  be finite and positive. The released coordinate remains the default.
- With omitted `sKnots`, asymmetric mode uses width-scaled negative anchors
  `[-w, -0.75w, -0.25w]` and the existing nonnegative anchors
  `[0, 0.25, 1.5, 3.5]`, with tiny endpoint padding for the strict fit support.
  Custom knots remain unchanged, and VIS rejects asymmetric mode.
- `backgroundmodel_BS(..., monotonic=True)` optionally represents each
  temporal component's spatial coefficients as a free baseline followed by
  nonnegative increments. The temporal spline, legacy coordinate, IRLS, and
  all row multipliers are unchanged.
- `monotonic_solver="normal"` reduces the weighted pixels to the small normal
  system before bounded fitting and is the operational default. `"full"`
  solves the equivalent bounded problem from all weighted rows as a reference.
  The ordinary released solver remains the default with `monotonic=False`.
- Monotonic BVLS retains `max_iter=1000` and accepts a result only when
  `result.success` is true. Failure raises by default. The explicit
  `monotonic_fallback=True` option warns, solves the failed IRLS iteration with
  the ordinary solver, and does not attempt another monotonic solve.
- Scalar `dgsolver_flag` is always saved: 0 for a standard requested fit, 1
  for monotonic success throughout, and 2 for monotonic failure followed by
  standard fallback. It is independent of optional diagnostics.
- `examples/scripts/validate_background_revision.py` accepts explicit native
  WIC, SI12, and SI13 orbit NetCDF files and runs the six cumulative Step 7
  variants. It writes representative four-field comparison figures and a CSV
  of fit time, robust summary metrics, iteration count, convergence, and the
  solver flag. Monotonic fallback is explicitly disabled so validation exposes
  solver failure. It is a comparison tool, not a validation result.
- Existing detector variance is preserved and checked on the fitted support.
  SI variance can be reconstructed when absent because SI is not reflattened.
  WIC variance cannot be reconstructed from the normally reflattened NetCDF
  image; WIC input without stored pre-reflattening `img_variance` is rejected
  with instructions to generate it through `read_idl`. Stored validation
  variance must identify `poisson_offset=0.5`; WIC must additionally carry the
  loader's `calculated_from="img before optional WIC reflattening"` provenance.
- `backgroundmodel_BS(..., diagnostics=True)` adds only the scalar
  `dgfit_iterations` and `dgfit_converged` fields. They are absent by default.

## Verification

The focused implementation and script-helper tests pass 59 tests. They cover the
dgsigma bin structure, known variance relation, and population invariance,
plus equal within-bin weights, equal bin influence, frame independence, and
excluded/empty-bin behavior, the legacy-coordinate scope, iteration-0 scale
retention, the released iteration-1 default, detector variance formulas,
calibration-table selection, loader opt-in behavior, masks and metadata, and
unsupported calibration states, plus constant-variance equivalence, reduced
high-variance influence, dgsigma separation, invalid fit variance, monotonic
coefficient ordering and profiles, normal/full equivalence, preserved ordinary
behavior, all solver flags, strict failure, warning/fallback behavior, one-time
monotonic attempts after failure, cumulative variant construction, frame
selection, plotting scales, stored WIC variance preservation and
rejection when absent or untraceable, SI variance attachment, and opt-in
scalar convergence diagnostics. CLI help passes.

The 46 focused BS tests also pass after adding asymmetric legacy coverage.
They verify the exact transform, unchanged legacy default, default and custom
width behavior, default and custom knots, invalid-width and VIS rejection,
and a combined fit using all revised BS options.

## Current boundary and next action

No geometry models beyond the single opt-in asymmetric legacy transform,
background covariance, diagnostic histories, or generated figures are present
on this branch. Complete-orbit execution is
blocked because `/home/bing/dtu_server` is not mounted, the available
icBuilder example IDL files are partial sequences, and a complete WIC orbit
must be regenerated with stored pre-reflattening variance. Package defaults
remain unchanged pending complete-orbit validation.
