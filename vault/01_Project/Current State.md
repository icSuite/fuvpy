# Current State

Last reviewed: 2026-09-12
Branch: `background_model_asymmetry`
Implementation status: BS-only directional model implemented; SH/SBS retired
Validation status: focused verification complete; complete-orbit comparison remains

## Production-readiness audit

The 2026-09-12 audit found no new numerical defect in the selected directional
BS path, but complete IMAGE processing must wait. The icBuilder orbit runner
still calls the old loader and obsolete SH stage, production BS damping has
not been selected, the selected model has not been checked on complete orbits,
and the runner is not interruption-safe. The loader also omits variance
provenance required by its own WIC validation path. Full findings are in
`vault/03_Validation/Production Readiness Audit - 2026-09-12.md`.

## BS-only architecture

`backgroundmodel_SH()` and `backgroundmodel_SBS()` have been removed together
with their dedicated helper, spherical-harmonic utility, accessor/export
paths, and active examples. The maintained background-removal path is now
`backgroundmodel_BS()` with optional directional asymmetry, and the corrected
image is `dgimg`. Historical NetCDF files containing SH fields remain readable,
but fuvpy no longer creates those fields.

## Historical comparison data

`examples/scripts/extract_legacy_background.py` now extracts only `dgmodel`,
`dgimg`, `dgweight`, and `dgsigma` at the bundled IDL example times from the
large historical orbit files in `~/Downloads/data`. The compressed references
are stored under `data/legacy_background/` and loaded by
`examples/scripts/compare_models.py` as the genuine historical result.

WIC, SI12, and SI13 each contain all 58 example frames. The three compressed
files occupy about 28 MB in total. The extraction script writes through a
temporary file and replaces each compact reference only after saving succeeds.

The historical BS fields were produced from complete northern-hemisphere
orbits with `x=cos(SZA)/cos(DZA)`, cubic spatial splines at
`[-3.5, -0.25, 0, 0.25, 1.5, 3.5]`, five quadratic temporal knots,
`tukeyVal=5`, `stop=0.01`, and absolute `dampingVal=1e-3`. There was no
monotonic constraint. Orbit 968 used 326 frames from 14:00 to 01:09; the
bundled example contains only three frames at 23:37--23:41. The saved result
can therefore be used as the genuine historical output, but it is not a
same-input comparison with models refitted to the short example sequence.

## Directional-model extraction

This branch starts from `background_model_revision` commit `1f6492f`. The
experimental `2D_model` branch is preserved separately at `1556ff1` and will
not be merged wholesale.

The only directional model selected for clean implementation is
`B0(x) + T(SZA)*S*Bs(x)`, where `x=cos(SZA)/cos(DZA)`, `S` is the pole-safe
solar-view sine factor, and `T` is optionally the Spencer taper. Both the sine
and cosine solar-view factors must be calculated and stored for visualization,
although the cosine factor will not enter the fit. The implementation plan is
in `vault/02_Algorithm/Directional Background Model Implementation Plan.md`.
The evidence and rejected alternatives from `2D_model` are summarized in
`vault/03_Validation/2D Model Experiments - Summary.md`.

The clean implementation now provides:

- `read_idl(..., viewing_geometry=True)`, which reconstructs DZA and stores
  both pole-safe solar-view factors;
- the packaged Spencer spherical-Chapman profile;
- `spatial_model="legacy_x_azimuth"` with configurable
  `directional_x_knots` and `spencer_taper`;
- separate `dgmodel_base` and `dgmodel_sine` outputs; and
- base-only monotonicity when `monotonic=True`.

The combined model reuses the revised IRLS, binned `dgsigma`, spatial and
measurement weighting, damping, convergence, and support masking. No
base-first fitting sequence or discarded 2D model machinery was introduced.

Twenty-one focused tests pass. The complete suite reports 100 passes and the
same three pre-existing failures documented on the parent branch. No complete
orbit has been run.

## Audit checkpoint

The 2026-09-09 read-only audit found that the core tensor-product spline design and monotonic transform are internally consistent, but three implementation defects must be repaired before scientific validation: output fields do not preserve actual fit support, unknown keywords are silently ignored, and asymmetric default knots are overwritten by a hard-coded experimental list.

The audit also corrected the interpretation of the recent Tukey change. Using `sqrt(Tukey)` as the row multiplier is conventional IRLS, but the Ohma paper and released implementation explicitly use the Tukey value itself inside a formulation with `W**2`. The new behavior is therefore a method revision, not an unambiguous bug fix. The same distinction applies to revised `1/sqrt(N)` bin-population weighting and dimensionless damping. Paper, later-package, and revised configurations must be made separately reproducible.

Full findings and the repair order are in `vault/03_Validation/BS Code Audit - 2026-09-09.md`. On 2026-09-10 the user clarified that improvement, not exact reproduction, is the objective. The conventional `sqrt(Tukey)` formulation remains the revised candidate. The production support-mask and keyword defects were fixed, the BS damping scale was made robust to unsupported columns, and `dgsigma` now uses the squared spatial row multiplier. Asymmetric-coordinate and residual-norm work are deferred.

## Confirmed implementation

- The released coordinate, knots, and noise-model choices remain the defaults.
  Robust row weighting and damping scaling have been revised on this branch,
  so the fitted result is no longer numerically identical to the release. The
  audit found that these are scientifically meaningful method changes, not
  unambiguous corrections to the published implementation.
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
  per-frame spatial sampling weights while the later-package Gaussian weighting
  remains the default. The paper-era implementation instead used `1/N` row
  multipliers. The option is limited to sensors using the legacy FUV
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
- WIC loading also accepts `wic_variance_model="correction_aware"`. This uses
  `10.65*a` instead of the constant 10.65 AD coefficient, where `a` is the
  reconstructed mission-gain and flat-field multiplier. The default remains
  `"reference"`, and the selected model is stored on `img_variance`. This
  scales conditional detector variance; it does not describe uncertainty in
  the calibration multiplier.
- `backgroundmodel_BS(..., input_variance="img_variance",
  variance_weighted=True)` applies `1/sqrt(variance)` only to the background
  inversion. Robust Tukey, spatial, and measurement row multipliers remain
  separate and visible. Empirical `dgsigma` and Tukey residual scaling do not
  use measurement variance.
- Measurement variance must be finite and strictly positive at every fitted
  pixel. BS scales zeroth-order damping as
  `dampingVal` times the median supported positive diagonal of the initial
  weighted normal matrix, freezing that scale through IRLS. A uniform change in row precision
  therefore no longer changes the relative damping strength. The numerical
  meaning of historical `dampingVal` settings has changed.
- The revised fit treats Tukey values as conventional objective weights and
  applies `sqrt(Tukey)` as the row multiplier. The Ohma paper and released code
  use the Tukey value directly in a formulation with squared row weights. Both
  definitions must remain reproducible until the revised choice is validated.
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

The revised Tukey and damping calculations pass 34 focused BS tests, including
ordinary and monotonic paths, measurement weighting, spatial weighting,
`dgsigma`, freeze behavior, and diagnostics. A constant measurement variance
with nonzero dimensionless damping reproduces the unweighted result, confirming
that the initial-`GTG` scale removes arbitrary global precision scaling.

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

At an earlier checkpoint, 46 focused BS tests passed after adding asymmetric
legacy coverage. The current full-suite audit supersedes that statement: 73
tests pass and three fail. One failure confirms that an experimental hard-coded
knot list now overrides the documented asymmetric knots; two failures concern
missing measurement-variance metadata.

Four focused tests for the reference and correction-aware WIC variance paths
pass. A broader set of 36 loader, BS measurement-weighting, and validation
tests passes when two pre-existing metadata tests are deselected. Those two
tests fail because the committed loader no longer writes the `excludes` and
`reflattening_assumption` attributes expected by the tests; this predates the
correction-aware option. A real 2000-08-28 WIC IDL file also loaded through
both paths; the correction-aware/reference variance ratio had 5th, 50th, and
95th percentiles of 1.02, 1.12, and 1.51.

## Current boundary and next action

Repair the variance metadata and update the production runner to use the
BS-only path, then run a 3--5-orbit complete BS pilot and select dimensionless
BS damping. Do not start complete IMAGE processing before that pilot is
accepted. Defer asymmetric-coordinate and residual-norm work.

No geometry models beyond the single opt-in asymmetric legacy transform,
background covariance, or diagnostic histories are present on this branch.
Complete-orbit validation is deliberately deferred until the audit defects and
reference configurations are repaired. Package defaults must not be changed
on the basis of the current comparisons.

Earlier SH diagnostics are retained in historical session and validation
notes, but no longer define active implementation work.
