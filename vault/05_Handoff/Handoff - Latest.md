# Handoff - Latest

Last updated: 2026-09-09
Branch: `background_model_revision`
Checkpoint: Steps 1 through 6 complete
Validation: script prepared; complete-orbit execution pending

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

The BS fit now has an explicit opt-in measurement weighting path:
`input_variance` names the Dataset variance variable and
`variance_weighted=True` applies `1/sqrt(variance)` to each fitted data row.
The calculation visibly combines robust, spatial-sampling, and measurement
row multipliers. Measurement variance is deliberately absent from both
empirical dgsigma regressions and Tukey residual scaling. Fitted variance must
be finite and strictly positive. The default remains unweighted.

Absolute measurement precision changes the balance between data and the
existing regularization term. Any nonzero `dampingVal` therefore requires a
new scale assessment rather than assuming its released tuning transfers.

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

## Next action

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
reusable combined-validation script. Scientific validation remains pending
because `/home/bing/dtu_server` is not mounted and the complete WIC input must
include correctly generated pre-reflattening variance; package defaults and
downstream production guidance are unchanged.
