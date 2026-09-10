# Handoff - Latest

Last updated: 2026-09-09
Branch: `background_model_revision`
Checkpoint: Steps 1 through 6 complete
Validation: read-only code audit complete; repairs required before complete-orbit execution

## Audit checkpoint

The 2026-09-09 BS audit is recorded in `vault/03_Validation/BS Code Audit - 2026-09-09.md`. It verified the tensor-product design and monotonic transform, but confirmed three defects: returned fields exceed the actual fit-support mask, unrecognized keywords are silently ignored, and `_asymmetric_legacy_knots()` overwrites its calculated knots with a hard-coded list. On 2026-09-10 the user deferred the asymmetric-coordinate defect because that path remains experimental.

The audit also established that the recent `sqrt(Tukey)` change is conventional IRLS but differs from the explicit `W**2` convention in Ohma et al. and the released implementation. Revised bin-population weighting and dimensionless damping likewise change the method. The validation script's `released` label is therefore not currently accurate.

On 2026-09-10 the support mask and keyword handling were fixed. Historical `dzalim` and `minlat` keywords are retained only as warned compatibility arguments for the Ohma paper scripts. Both BS and SH now scale damping from numerically supported positive `GTG` diagonal entries, and the empirical BS `dgsigma` calculation uses `Tukey * spatial_multiplier**2`. Measurement variance remains excluded from `dgsigma` intentionally.

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

WIC loading additionally supports
`wic_variance_model="correction_aware"`. It uses `kappa_i=10.65*a_i`, where
`a_i` is the reconstructed mission-gain times flat-field multiplier, in
`variance=kappa_i*img + poisson_offset*kappa_i**2`. The constant 10.65
coefficient remains the `"reference"` default. The selected model is stored
on `img_variance`. This remains conditional on the multiplier being correct
and does not add flat-field or gain uncertainty.

The BS fit now has an explicit opt-in measurement weighting path:
`input_variance` names the Dataset variance variable and
`variance_weighted=True` applies `1/sqrt(variance)` to each fitted data row.
The calculation visibly combines robust, spatial-sampling, and measurement
row multipliers. Measurement variance is deliberately absent from both
empirical dgsigma regressions and Tukey residual scaling. Fitted variance must
be finite and strictly positive. The default remains unweighted.

BS and SH now use dimensionless zeroth-order damping:
`R=dampingVal*scale*I`, where `scale` is the median supported positive diagonal
of the initial weighted normal matrix. It remains fixed while IRLS updates the Tukey weights. A
uniform change in row precision therefore does not alter relative damping,
but historical numeric `dampingVal` settings do not transfer directly.

The stored Tukey fields remain conventional objective weights. Both fits use
their square roots as row multipliers, so the least-squares objective contains
the Tukey weight once rather than squaring it again. Spatial and measurement
terms remain direct row multipliers.

The optional `asymmetric_legacy=True` coordinate starts from
`x=cos(SZA)/cos(DZA)`, preserves `x>=0`, and maps negative values to
`nightside_width*tanh(x/nightside_width)`. The width defaults to 0.2 and must
be finite and positive. When knots are omitted, the negative anchors are
`[-w, -0.75w, -0.25w]` and the remaining anchors stay
`[0, 0.25, 1.5, 3.5]`; only tiny endpoint padding is added for strict support.
Custom knots are respected. The released coordinate remains the default, and
VIS rejects asymmetric mode. The transform uses the existing BS design,
noise, spatial weighting, robust weighting, measurement weighting, and
monotonic machinery without alternative paths.

The optional monotonic BS fit uses a free first spatial coefficient plus
nonnegative increments for every temporal basis component. This matches the
existing spatial-outer, temporal-inner coefficient ordering and leaves the
temporal B-spline unchanged. `monotonic_solver="normal"` first reduces the
weighted image rows to the coefficient-sized normal system;
`monotonic_solver="full"` is the full-row bounded reference. Both retain the
same robust, spatial, and measurement row multipliers. The default remains the
ordinary unconstrained solver. BVLS keeps `max_iter=1000`, and `result.success`
is the criterion for accepting the bounded result.

Monotonic failures raise `RuntimeError` by default. With the explicit
`monotonic_fallback=True` option, failure instead emits a `RuntimeWarning` and
the ordinary solver handles that IRLS iteration and all remaining iterations;
the failed bounded solve is not retried. Every fit saves scalar
`dgsolver_flag`: 0 standard requested, 1 monotonic succeeded throughout, and 2
monotonic failed and standard fallback was used. This output does not depend
on diagnostics.

`examples/scripts/validate_background_revision.py` is a concise Step 7
comparison tool. It requires explicit complete native-grid WIC, SI12, and SI13
orbit NetCDF paths and runs the six cumulative variants from released through
monotonic. Existing `img_variance` is preserved and checked on fitted pixels.
Missing SI variance is attached with `poisson_offset=0.5`; missing WIC variance
is rejected because the normally stored WIC image is already reflattened and
cannot reproduce the loader's pre-reflattening calculation. A suitable WIC
orbit must be loaded and saved with
`read_idl(..., measurement_variance=True, poisson_offset=0.5)`.
Stored validation variance must identify `poisson_offset=0.5`. WIC variance
must also retain the loader attribute
`calculated_from="img before optional WIC reflattening"`; otherwise the script
rejects it rather than accepting variance with ambiguous provenance.

The script writes representative figures showing `dgmodel`, `dgimg`,
`dgsigma`, and `dgweight`, plus a CSV containing fit time, robust summary
metrics, IRLS iteration count, convergence, and `dgsolver_flag`. Fallback is
explicitly disabled in the monotonic validation variant so a bounded-solver
failure remains visible. Iteration count and convergence are supplied by the
minimal opt-in `backgroundmodel_BS(..., diagnostics=True)` fields
`dgfit_iterations` and `dgfit_converged`; neither field is returned by default.
No diagnostic histories were added. The script does not run the SH model or
change package defaults.

## Verification

The audit result supersedes older passing-test statements below where they
conflict: after the 2026-09-10 repairs, the current complete suite has 79 passes and three failures. The BS
failure is the hard-coded asymmetric-knot override; the other two concern
missing measurement-variance metadata.

The revised Tukey and dimensionless-damping calculations pass 34 focused BS
tests covering ordinary and monotonic solves, measurement and spatial weights,
`dgsigma`, freeze behavior, and diagnostics. The new nonzero-damping test
confirms invariance to a uniform measurement-precision scale. A synthetic
BS-to-SH smoke test completed with finite SH output. A broader 76-test run had
73 passes and three failures unrelated to these edits: two known missing
measurement-variance metadata fields and the current experimental hard-coded
asymmetric-knot override.

`pytest -q tests/test_background_dgsigma.py
tests/test_background_spatial_weighting.py
tests/test_background_dgsigma_freeze.py
tests/test_measurement_uncertainty.py
tests/test_background_measurement_weighting.py
tests/test_background_monotonic.py
tests/test_background_revision_validation.py
tests/test_background_diagnostics.py` passes 59 tests. In addition to Steps
1--3, tests cover SI and WIC equations and offsets, representative mission and
flat-field selection, dimensions, masks, metadata, rejected calibration
states, pre-reflattening calculation, and unchanged default loader behavior.
Step 5 tests confirm zero-damping equivalence for constant variance, reduced
influence of a high-variance observation, separation from dgsigma estimation,
and rejection of invalid fitted variance.
Step 6 tests cover coefficient ordering, monotonic legacy-coordinate profiles,
normal/full equivalence with combined row weights and damping, unchanged
ordinary behavior, all three solver flags, strict solver failure, explicit
warning/fallback behavior, and no repeated monotonic attempt after fallback.
The script tests cover the exact cumulative variant sequence, representative
frames, shared robust limits, stored WIC variance preservation, explicit WIC
rejection when variance is absent, invalid, or lacks the required offset and
pre-reflattening provenance, SI attachment, and required CLI arguments.
Diagnostics tests cover default absence and opt-in scalar converged values.
`python examples/scripts/validate_background_revision.py --help` also passes.
No orbit model comparison was run.

The WIC variance option has four passing focused tests. A broader group of 36
loader, BS measurement-weighting, and validation-script tests passes after
deselecting two pre-existing metadata tests. Those tests fail because the
current committed loader no longer writes `excludes` and
`reflattening_assumption`; this was not changed as part of the WIC option. A
real 2000-08-28 WIC file also loaded successfully with both models; the
correction-aware/reference variance ratio had 5th, 50th, and 95th percentiles
of 1.02, 1.12, and 1.51.

A focused compatibility check on the bundled 15-frame WIC sequence ran the
revised variance-weighted BS fit, the existing SH fit, `Dataset.to_netcdf()`,
and `xarray.load_dataset()` successfully. `img_variance`, all revised BS
fields, `dgsolver_flag`, and the SH products survived the round trip. The SH
stage remains scientifically unchanged: it normalizes with `dgsigma` and
inherits `dgweight`, but does not apply `img_variance` as measurement
precision. The result therefore confirms structural and serialization
compatibility, not a fully variance-aware SH inversion.

The first SH comparison script confirmed that the apparently unchanged image
is primarily a display-scale issue: `shimg` equals `dgimg - shmodel` exactly,
but the smooth SH correction is much smaller than the bright auroral range.
For WIC at 2000-05-16 03:03:32, the SH correction was about 9.7% of the
`dgimg` RMS and approximately -100 to +142 counts. The script's ordinary and
"variance" SH rows are identical because `backgroundmodel_SH` silently
ignores `variance=True` and other unknown BS-only keywords. A future SH API
revision should validate keywords rather than silently accepting a comparison
that did not occur.

The August 28 WIC sequence exposed stronger SH shrinkage than the display
alone suggested. Median `dgsigma` was about 341 counts and the equal-area bins
held a median 129 pixels. SH used `1/N` as a row multiplier, so the normal
equations applied `1/N^2` per pixel. The resulting data normal matrix had
median diagonal `1.34e-5`, below the fixed `1e-4` damping, and only 5 of 90
data eigenvalues exceeded the damping. The row multiplier is now corrected to
`1/sqrt(N)`. BS and SH now both scale damping from their initial `GTG` and use
`sqrt(Tukey)` as the robust row multiplier. The earlier fixed-damping sweep is
therefore obsolete; complete-orbit comparison must select a dimensionless
damping value.

`pytest -q tests/test_background_asymmetric_legacy.py
tests/test_background_dgsigma.py tests/test_background_spatial_weighting.py
tests/test_background_dgsigma_freeze.py
tests/test_background_measurement_weighting.py
tests/test_background_monotonic.py tests/test_background_diagnostics.py`
passes 46 tests. Asymmetric tests cover the exact mapping, released-coordinate
identity, default width and knot anchors, a custom width, invalid widths, VIS
rejection, preserved custom knots, and a combined revised-option fit.

## Next action

Run the complete-orbit comparison and select dimensionless damping. Keep `sqrt(Tukey)` as the revised candidate and make comparison labels truthful; exact paper reproduction is not the development objective. Defer asymmetric-coordinate and residual-norm work.

Mount the complete native-grid orbit data, generate and save the complete WIC
orbit with pre-reflattening `img_variance`, then run the script with explicit
paths, for example:

```bash
python examples/scripts/validate_background_revision.py \
  --wic /path/to/complete_wic_orbit.nc \
  --si12 /path/to/complete_si12_orbit.nc \
  --si13 /path/to/complete_si13_orbit.nc \
  --output-dir /tmp/fuvpy_background_revision_validation
```

Do not use the partial icBuilder example IDLs for scientific validation and do
not change package defaults before the complete-orbit comparison is accepted.

## Portfolio impact

The clean background-model revision has completed Steps 1 through 6 and has a
reusable combined-validation script plus one isolated asymmetric legacy
coordinate option. Scientific validation remains pending
because `/home/bing/dtu_server` is not mounted and the complete WIC input must
include correctly generated pre-reflattening variance; package defaults and
downstream production guidance are unchanged.
