# Decision Log

## 2026-09-12 — Remove obsolete examples and rejected model paths

**Decision:** Remove the article-specific examples and generated comparison
figures, and retain only the three scripts and detector sequences used for
current BS validation. Also remove the rejected compressed-nightside
`asymmetric_legacy` option. The only retained asymmetry is
`spatial_model="legacy_x_azimuth"`.

**Rationale:** The article examples describe an algorithm that is no longer
the active method. Keeping generated outputs and two meanings of “asymmetry”
made the repository larger and the scientific interface less clear.

## 2026-09-12 — Retire the SH and SBS background paths

**Decision:** The active fuvpy background-removal architecture consists only
of `backgroundmodel_BS()`, with optional directional asymmetry. Remove
`backgroundmodel_SH()`, `backgroundmodel_SBS()`, their dedicated helper and
utility code, and their public accessor/export paths. New corrected-image
workflows use `dgimg`; historical files containing `shimg` remain readable as
ordinary xarray datasets.

**Rationale:** Neither the secondary spherical-harmonic correction nor the
alternative SBS model belongs to the selected production architecture.
Removing them makes the scientific processing path explicit and prevents
obsolete stages from being invoked accidentally.

## 2026-09-11 - Extract one directional model from the experimental branch

**Decision:** Preserve `2D_model` at `1556ff1` as the experimental record and
start `background_model_asymmetry` from `background_model_revision` at
`1f6492f`. Reimplement only `B0(x) + T(SZA)*S*Bs(x)`. Calculate and save both
solar-view sine and cosine factors, although only the sine factor enters the
fit. Expose the Spencer taper as an option.

**Rationale:** The experimental branch contains many rejected model variants
and supporting switches. Selective reimplementation retains the most
promising directional correction without turning the revised production code
into a generic model framework. The unused cosine factor remains valuable for
visual inspection of the viewing geometry.

## 2026-09-03 — Attach detector variance during opt-in IDL loading

**Decision:** Package unchanged copies of the FUVVIEW3 mission and flat-field
tables and let `read_idl(..., measurement_variance=True)` reconstruct the
combined calibration multiplier and add `img_variance`. Preserve ordinary
loader behaviour by default. Use the WIC reference law by default within the
opt-in path and require an explicit `wic_variance_model="correction_aware"`
choice for the additional multiplier.

**Rationale:** Detector variance must be evaluated from the calibrated input
image before background subtraction. Keeping the multiplier and variance in
the loaded dataset makes the calculation inspectable and lets the background
fit consume `input_variance="img_variance"` without duplicating sensor physics.
The explicit WIC choice avoids presenting the unresolved stellar-sensitivity
interpretation as settled.

Last reviewed: 2026-09-03

## 2026-09-03 — Treat WIC `IMAGE` as flag-1-equivalent corrected AD counts

**Decision:** Use the working interpretation that the saved WIC `IMAGE` field
contains the corrected-AD product produced by the flag-1 calibration path.
Treat flag 3 as the label for the additional `CIMAGE` background-subtracted
product, not as evidence that `IMAGE` itself has been converted to Rayleigh.

**Rationale:** The recovered FUVVIEW3 source defines flag 1 as corrected counts
and flag 2 as Rayleigh. WIC files contain both a nonnegative `IMAGE` field and a
signed `CIMAGE` field, while fuvpy converts `IMAGE` to kR using the published
612.6 AD/kR sensitivity. This remains an explicit working assumption because
the producer-side flag-3 branch is not present in the recovered source.

## 2026-09-02 — Stop the background-only variance model at Step 3

**Decision:** Do not propagate the Step 2 background-only effective variance
laws into downstream ensembles or icBuilder. The agreed statistical gate
failed for WIC and SI13. Treat a model with SZA, mission-time, or orbit-specific
dependence as a new model-development cycle requiring explicit approval and
new held-out validation, not as an automatic continuation to Step 4.

**Rationale:** Four unused complete-orbit fits gave pooled standardized spreads
of 0.611 for WIC, 1.010 for SI12, and 0.645 for SI13. WIC and SI13 also retained
large background, SZA, and orbit dependence. SI12 alone passed, which is
insufficient for the joint three-sensor processing chain.

## 2026-09-02 — Start uncertainty estimation at the saved IDL product

**Decision:** Do not seek Level-1 UDF files or modify FUVVIEW3 to export
additional detector statistics. Define the required input as an effective
observation variance in the corrected-count space of the saved IDL images.
Recovered FUVVIEW3 calibration tables and instrument literature may constrain
or check this estimate, but the IDL files are the earliest data product in
scope.

**Rationale:** Reconstructing the original detector pipeline is not feasible
for this project. SI has a physically useful photon-counting interpretation,
but WIC's intensified-CCD response and the undocumented calibration flag in
the saved WIC files prevent a uniform exact inversion to raw detector noise.
The common feasible route is therefore a robust, sensor-specific estimate from
background-dominated post-FUVVIEW residuals, explicitly labelled as effective
observation variance rather than pure detector variance.

## 2026-09-02 — Condition on the selected background model

**Decision:** For this uncertainty component, assume that the selected fuvpy
background model is correct. Quantify observation-noise propagation and fit
uncertainty conditional on that model; do not attempt to estimate
background-model discrepancy.

**Rationale:** Model discrepancy is not identifiable from the current fit
covariance and would require a separate comparison of background model
families. It should not block the narrower, explicitly conditional product.
The immediate blocker is instead the variance of the detector image supplied
to the fit.

## 2026-09-02 — Defer robust covariance and downstream integration

**Decision:** Keep the fixed-weight covariance as an experimental conditional
spread/correlation product, and retain complete repeated fits as the WIC
reference. Do not add a sandwich estimator or integrate the result into
icBuilder until a calibrated detector observation-variance model exists.

**Rationale:** Complete-orbit tests show that the fixed-weight calculation
captures marginal spread and short-lag correlation reasonably under stated
perturbations. The main discrepancy is a nonlinear repeated-fit mean shift,
which a covariance or sandwich correction does not itself encode. The current
`dgsigma**2` and count-scale diagnostic are not calibrated detector variances.

## 2026-09-02 — Keep analytic and repeated-fit products distinct

**Decision:** Expose fixed-weight coefficient covariance as an experimental,
conditional sampling approximation. Use complete repeated fits as the
reference when robust weights or the dependent WIC BS-to-SH chain matter. Do
not construct WIC uncertainty by adding separate BS and SH variances.

**Rationale:** The analytic expression is correct only after freezing the fit
mask, weights, basis, and regularization. Repeating the complete chain retains
data-dependent weights, the early-fit spread calculation, and BS/SH
dependence. Complete-orbit evidence confirms that the two products answer
different questions and should remain distinct.

## 2026-09-01 — Isolate background-uncertainty development

**Decision:** Develop fitted-background uncertainty on the local
`background_uncertainty` branch. Preserve the existing background estimates
and default API while the uncertainty definition is tested.

**Rationale:** icBuilder needs observation uncertainty before coregistration
and nonlinear propagation, but changing the selected background model at the
same time would confound a mean-processing change with an uncertainty change.

## 2026-09-01 — Keep uncertainty sources distinct

**Decision:** Treat residual/count noise, conditional fitted-background
uncertainty, and background-model discrepancy as separate quantities. Do not
label `dgsigma` or the inverse normal matrix as total posterior uncertainty.

**Rationale:** The robust fit deliberately leaves auroral signal in the
residuals, its weights are data-dependent, and regularization is not yet a
calibrated probability prior. Systematic model and calibration errors are not
represented by coefficient covariance.

## 2026-09-04 — Keep measurement-weighted BS fitting explicit and opt-in

**Decision:** Add `variance_weighted=True` to the BS model. Require an explicit
input variance, use (1/\sqrt{v_i}) in the least-squares rows, and calculate
Tukey weights from variance-standardized residuals. Preserve the original fit
as the default and reject non-positive variance at fitted pixels.

**Rationale:** Merely propagating variance after the fit cannot change the
background estimate or robust weights. The requested heteroscedastic fit must
use measurement variance inside IRLS. Keeping it opt-in separates this
scientific change from established fuvpy behavior. Its interaction with the
heuristic damping value still requires empirical assessment.

**Refinement:** Preserve measurement-standardized Tukey weighting as an option,
but use `tukey_scale="empirical"` for the first WIC comparison. Coherent
dayside model residuals can exceed several detector-noise standard deviations
without being isolated bad pixels. Also expose `tukey_weights=False` so robust
rejection can be tested independently from measurement weighting.

**Low-count rule:** Extend `poisson_offset` to WIC equivalent output quanta.
Use `poisson_offset=0.5` for variance-weighted experiments, giving the reference
model a floor of \(0.5(10.65)^2\) AD squared. Retain zero as the API default so
existing loading behavior does not change.

**Geometry test:** Expose an optional `sza_max` in the BS fit so nightside SZA
can be clipped experimentally before calculating `cos(SZA)/cos(DZA)`. Preserve
unclipped behavior by default and never overwrite the dataset's SZA values.

**Noise-model correction:** Add `noise_model="binned"` using one paired
residual-variance and mean-model-squared value per populated geometry bin.
Retain `noise_model="legacy"` as the default until full-orbit comparison is
complete. This corrects pixel-level duplication of bin statistics without
changing the published bin layout, nonnegative coefficients, or early scale
freezing.

**Experimental controls:** Expose `gaussian_weights` and
`dgsigma_freeze_iteration` rather than silently replacing the published/released
fit variants. Keep the current defaults (`True` and 1) for reproducibility;
use `False` and 0 to test an unweighted spatial fit with the residual scale
frozen after its initial estimate.

## 2026-09-10 — Prioritize the revised model over exact reproduction

**Decision:** The purpose of `background_model_revision` is to improve the BS
background estimate used by the conductance pipeline, not to reproduce the
Ohma implementation exactly. Retain the conventional `sqrt(Tukey)` row
multiplier as the revised candidate. Keep the published convention identifiable
only where needed to make comparisons interpretable.

**Scope:** Treat `asymmetric_legacy` as experimental and defer its hard-coded
knot defect. Defer correction of the stored residual/model norms until those
diagnostics are needed for damping selection or an L-curve. The immediate
production issues are the output support mask, silent keyword handling, and a
safe nonzero damping scale.

## 2026-09-13 — Fix dgsigma after the initial fit

**Decision:** `backgroundmodel_BS()` always estimates `dgsigma` after IRLS
iteration 0 and then holds it fixed while the Tukey weights and model continue
to update. Remove the `dgsigma_freeze_iteration` option.

**Rationale:** Iteration 0 is the intended revised calculation and no
production use remains for the experimental freeze-iteration variants. A
single behavior makes the fit and its metadata easier to interpret.

## 2026-09-13 — Fix the noise model and spatial weighting

**Decision:** `backgroundmodel_BS()` always uses the binned `dgsigma` model and
bin-population spatial weighting. Remove the `noise_model` and
`spatial_weighting` options, the pixel-expanded legacy noise model, and the
Gaussian weighting path.

**Rationale:** These are retained corrections rather than unresolved model
choices. Removing the historical alternatives narrows the production function
and prevents accidental processing with settings already rejected for the new
archive.

## 2026-09-13 — Remove unused BS interface controls

**Decision:** `backgroundmodel_BS()` always returns a copy, always stores
convergence diagnostics, uses `img_variance` whenever variance weighting is
enabled, and uses the low-cost normal-system monotonic solver in production.
Remove the public `inplace`, `diagnostics`, `input_variance`, and
`monotonic_solver` selectors. Also reject the obsolete `dzalim` and `minlat`
fit arguments; geometry filtering belongs to `read_idl()`.

**Rationale:** None of these were unresolved scientific choices. Fixing their
behavior makes production calls shorter and prevents configurations that were
ignored, redundant, or retained only for superseded comparisons. The full
monotonic solver remains internal for numerical equivalence tests.
