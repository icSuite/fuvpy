# B-Spline Background Model Code Audit

Last reviewed: 2026-09-09
Branch: `background_model_revision`
Commit inspected: `afcaa7d` plus the tracked working-tree edits present on 2026-09-09
Scope: `backgroundmodel_BS`, its helpers, focused tests, validation scripts, relevant Git history, and the method described by Ohma et al. (2024)

## Executive conclusion

The tensor-product B-spline construction and the optional monotonic coefficient transform are internally consistent. The code nevertheless has three confirmed implementation defects that should be fixed before scientific validation:

1. returned fields do not preserve the actual fit-support mask;
2. unknown and obsolete keyword arguments are silently ignored; and
3. the asymmetric-coordinate knot helper discards its calculated knots and returns a hard-coded experimental list.

There is also a reproducibility problem. The current branch's nominal `released` comparison is not the published/released algorithm: Tukey weighting and damping have been changed globally, while the default Gaussian spatial weighting was itself introduced after the paper implementation as a test. These alternatives may be reasonable, but they must be named as distinct model definitions and compared explicitly.

No BS source code was changed during this audit.

## Post-audit implementation — 2026-09-10

The production-relevant support-mask and keyword defects are now fixed. All BS products are masked to the pixels that actually entered the fit. Unknown keywords raise immediately. The historical `dzalim` and `minlat` BS arguments were initially retained with a warning, then removed on 2026-09-13 after the obsolete paper examples were deleted; geometry limits now belong solely to loading.

The damping scale now uses the median of numerically supported positive entries of `diag(GTG)` in both BS and SH and raises if no coefficient is supported. The `dgsigma` regressions now receive `Tukey * spatial_multiplier**2`, matching the Tukey/spatial part of the revised inversion objective while continuing to exclude measurement variance deliberately.

## Confirmed high-priority defects

### 1. Returned fields exceed the fitted support

The fit mask requires finite image data and a spatial coordinate strictly inside the spline support. The output mask checks only `sza >= 0`.

Consequences:

- `dgmodel` can contain unconstrained spline extrapolation outside the knot range;
- `dgweight` can be 1 for a pixel that was never fitted;
- `dgmodel` and `dgweight` can remain finite where the input image was NaN;
- `dgsigma` is NaN on the same pixels, leaving internally inconsistent products.

A synthetic pixel with spatial coordinate 107.68, far outside the default upper knot at 3.5, returned `dgmodel=-31,804,645`, `dgweight=1`, and `dgsigma=NaN`. A finite-SZA pixel with `img=NaN` returned a finite model and weight 1.

The implementation should retain the three-dimensional fit-support mask. If predicting a background where the observation is missing is desired, prediction support must be defined separately from fit support; a fitted-data weight must never be reported there.

### 2. Unknown keywords are silently ignored

Recognized options are removed from `kwargs`, but the function never checks what remains. The calls below therefore produce identical results without warning:

```python
backgroundmodel_BS(images)
backgroundmodel_BS(images, dzalim=0, minlat=89, nonsense_option="ignored")
```

This affects reproducibility directly. The paper figure script still passes `dzalim` and `minlat`, although those options were removed from `backgroundmodel_BS`. The intended DZA exclusion now occurs indirectly because `read_idl(dzalim=...)` masks the coordinate fields, and icBuilder currently loads with `dzalim=75`. A dataset produced another way receives no equivalent protection.

Reject all unconsumed keywords. Separately decide whether geometry limits belong only in the loader or should be explicit BS options.

### 3. Asymmetric knots are overwritten

`_asymmetric_legacy_knots()` calculates width-dependent knots and then replaces them with a fixed six-element list. Changing `nightside_width` therefore changes the coordinate but not the knots. The fixed list also omits the zero knot and includes an unexplained positive knot at 0.05.

The corresponding focused test fails. Delete the hard-coded override and retain the documented width-dependent result.

## Reproducibility and method-definition issues

### Published Tukey convention versus conventional IRLS

Ohma et al. write the fit using `W^2` and define each Tukey value as `[1-u^2]^2`. The published/released code uses that value directly as a row multiplier, so its least-squares objective contains the fourth power of the bisquare expression.

The current working tree instead uses `sqrt(Tukey)` as the row multiplier. This is the conventional IRLS interpretation: the stored Tukey value enters the objective once. It is mathematically defensible, but it is a change from the published implementation, not an unambiguous coding correction.

Both definitions should remain reproducible during validation and be labeled clearly.

### Spatial data-coverage weighting

The paper and its implementation use `1/N` as a row multiplier for a geometry bin containing `N` samples. Because the normal equations square row multipliers, the total objective influence of that bin is proportional to `1/N`.

The revised `bin_population` option uses `1/sqrt(N)` and normalizes within each frame. This gives every populated bin equal total squared influence. That matches the stated equal-bin objective better, but again differs from the published implementation.

The current default `gaussian` path is also not the paper implementation. Git history shows it was introduced under the comment `TEST WITH NORMAL DIST` and then became the package default. At the edge of the default fit domain its row multiplier is about 15, corresponding to an objective multiplier above 200.

### Damping definition

The branch scales damping as

```text
dampingVal * median(diag(GTG)) * I
```

using the initial weighted normal matrix. This makes damping invariant to a uniform precision rescaling, which is useful, but historical values of `dampingVal` and the paper's L-curve choices no longer have the same numerical meaning.

The scale can also collapse to zero when unsupported coefficient columns make at least half of `diag(GTG)` zero. In a one-frame synthetic fit with `n_tKnots=5`, `dampingVal=1000` produced exactly the same model as `dampingVal=0`. The scale should use supported positive diagonal entries, and incompatible temporal configurations should be rejected.

### Reference comparisons are mislabeled

`validate_background_revision.py` calls an empty settings dictionary `released`, but global changes to Tukey and damping mean that this no longer reproduces the released fit. In `compare_bs_models.py`, every row called `Original` already uses bin-population weighting and freezes `dgsigma` at iteration zero.

A future validation should define explicit configurations such as:

- `paper`: paper-era `1/N` spatial row weights, paper-era Tukey row convention, absolute damping;
- `later_package`: Gaussian spatial weighting and the later package defaults;
- `revised`: binned `dgsigma`, equal-bin influence, chosen Tukey convention, measurement precision, dimensionless damping, and optional monotonicity.

## Medium-priority numerical and diagnostic issues

### The ordinary solver uses normal equations

The code solves `(G.T @ G + R)m = G.T @ d`. Forming `G.T @ G` squares the condition number and is less stable than solving the weighted, regularized least-squares problem as an augmented system with QR/SVD. The small synthetic audit case produced the same prediction, so this is a numerical risk rather than a demonstrated production failure.

### `dgsigma` weights did not match the fitted objective

At audit time, the revised fit objective was weighted by `Tukey * spatial_multiplier**2 * measurement_precision**2`, whereas the empirical `dgsigma` calculation received `Tukey * spatial_multiplier`. This was corrected on 2026-09-10 to use `Tukey * spatial_multiplier**2`. Measurement variance remains deliberately excluded because `dgsigma` is an empirical residual-spread model rather than detector measurement variance.

### Stored residual norms do not describe the fitted objective

`dgnorm_residual` uses only the final Tukey values. It omits spatial and measurement weights and therefore cannot serve as the residual norm for an L-curve of the revised inversion. `dgnorm_model` is an RMS coefficient magnitude rather than the conventional coefficient 2-norm. Existing L-curve code should not be used for the revised model until the norm definitions are aligned with its objective.

### Final model and final weights are offset by one IRLS update

Each iteration solves with the current weights and then calculates new weights. The returned model was therefore fitted with the preceding weights, while `dgweight` contains the newly calculated values. On the converged synthetic test this was negligible: one additional solve changed the model by `6.2e-7` relative at the default stopping tolerance. It can matter when the loop reaches the iteration limit or stops abnormally. Non-convergence should be an explicit warning/flag, and returned diagnostics should correspond to one internally consistent final state.

### Input validation is incomplete

Invalid values can fail late or return questionable output:

- zero `tukeyVal` causes division by zero and a later empty-fit error;
- negative damping can run to the iteration cap and still return a model;
- duplicate or unordered dates and incompatible knot/order choices yield cryptic SciPy errors;
- an unsupported sensor id produces an indirect exception;
- integer image arrays fail when NaN is assigned to arrays created with `ones_like` or `full_like`.

Validate the public scientific inputs before building the design matrix and create floating-point work arrays explicitly.

## Verified parts

- Spatial/temporal tensor-product coefficient ordering matches the reconstruction ordering.
- The monotonic transform correctly gives each temporal basis component a free baseline and nonnegative spatial increments.
- The normal and full bounded monotonic solvers agree in the focused weighted test.
- `1/sqrt(img_variance)` is the correct row multiplier for measurement precision.
- The revised `1/sqrt(N)` bin-population multiplier gives equal total squared influence to populated bins within a frame.
- The binned `dgsigma` regression removes the released code's accidental repetition of each bin statistic once per pixel and is closer to the method described in the paper text.
- With current icBuilder loading, the DZA cutoff is applied by masking coordinate fields before the BS fit. This does not remove the API problem for other callers.

## Test status

After the post-audit fixes, the complete local test suite produced 79 passes and three failures:

1. the confirmed asymmetric-knot overwrite;
2. missing `excludes` metadata on `img_variance`; and
3. missing `reflattening_assumption` metadata on WIC variance.

The last two failures concern loader metadata rather than the BS inversion, but they mean the branch is not currently test-clean.

## Recommended repair order

User disposition on 2026-09-10: exact reproduction is not the objective, the conventional `sqrt(Tukey)` formulation remains the revised candidate, asymmetric-coordinate work is experimental and deferred, and residual-norm correction is deferred until the diagnostic is needed.

1. ~~Fix the output support mask and reject unknown keywords.~~ Completed 2026-09-10.
2. Add validation for sensor, dates, spline settings, damping, Tukey scale, nonempty support, and positive supported `GTG` scale.
3. Make comparison labels truthful without making reproduction the development target.
4. Validate and lock the revised Tukey and spatial-weight conventions with objective-level tests.
5. ~~Align the Tukey/spatial `dgsigma` weights with the revised objective.~~ Completed 2026-09-10; measurement variance remains intentionally separate.
6. Defer the asymmetric knot, residual-norm, and normal-equation changes unless they enter the selected production path or validation.
7. Run the complete-orbit comparison and select damping.
