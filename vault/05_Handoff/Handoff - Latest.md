# Handoff - Latest

Last updated: 2026-09-08
Repository snapshot: `background_uncertainty` at `763de63`
Worktree state: uncertainty implementation, tests, documentation, and project
memory are uncommitted

## Current checkpoint

An experimental collapsed geometry is available through
`geometry_model="effective_illumination"`. It uses
`log(smooth_solar(SZA)) + viewing_exponent*log(van_Rhijn(DZA))`, so different
DZA values remain distinguishable at the terminator while the BS model stays
one-dimensional. Defaults are 130 km emission height, exponent 1, and a 5
degree twilight width; all are exposed. The legacy coordinate remains the
default. Effective illumination uses neutral spatial weighting unless the
caller explicitly selects bin-population weighting; the legacy Gaussian
weight is rejected because its scale is tied to the old coordinate.
The calculated detector-grid coordinate is now retained as
`dayglow_coordinate`, with its model parameters stored as attributes.
The coordinate and its default knots are translated so SZA 90 degrees and DZA
0 degrees define zero. The unit-effective-illumination knot therefore moves
to about 2.81. This changes labels but not the fitted model.
The additional `geometry_model="bounded_effective_illumination"` uses the
same physical terms with a shifted `log1p` transform. Its deep-night limit is
`-log(2)`, the nadir terminator remains zero, and the existing raw-log option
is preserved. `compare_bs_freeze_new.py` now exercises this bounded option
and writes to `figures/IRLS_iterations_bounded/`; all 129 example-frame
figures completed. Visual inspection of a representative WIC frame showed
some changes in the fit and Tukey masks, but no clear improvement in the
corrected image.
The further `geometry_model="linear_effective_illumination"` removes the
logarithmic compression while retaining the smooth solar term and finite
viewing factor. It normalizes the nadir terminator to zero and the nadir
subsolar point to one, with extra default knots across the dayside. The
comparison script preserves the user's manual raw-log knots, selects defaults
for the linear option, and writes separately to
`figures/IRLS_iterations_linear/`. All 129 example-frame figures completed;
the representative WIC correction remains broadly similar to the other
coordinate choices.
`geometry_model="asymmetric_legacy"` now provides the complementary test:
the dayside is exactly the legacy `cos(SZA)/cos(DZA)` coordinate, while the
negative side smoothly saturates through
`nightside_width*tanh(x/nightside_width)`. The exposed width defaults to 0.2.
Legacy Gaussian weighting remains available for this legacy-derived model,
and the original `legacy` implementation is untouched. Five focused tests
passed, including stored-background reconstruction, and all 129 comparison
figures completed under `figures/IRLS_iterations_asymmetric_legacy/`.
`compare_bs_freeze_new.py` now assigns `bin_population` to the correct
`spatial_weighting` setting instead of accidentally overwriting
`twilight_width`. Public geometry parameters also reject non-numeric values
with a direct error message.

For a fairer `dgsigma` comparison, `effective_illumination` and
`linear_effective_illumination` now use the legacy bin structure: one bin for
all negative coordinate values and 14 equal-width bins over the positive
range. This does not change their spline knots or spatial population weights.
Three focused bin/noise tests passed.

A separate `geometry_model="cos_sza"` now supplies the simplest comparison:
the coordinate is exactly `cos(SZA)`, with no DZA dependence or geometry
tuning parameters. It uses one broad nightside spline interval, one nightside
`dgsigma` bin, and 14 equal dayside `dgsigma` bins. Four focused tests passed;
all 129 comparison figures are in `figures/IRLS_iterations_cos_sza/`.

The separate `geometry_model="powered_legacy"` provides a controlled bridge
between the pure-SZA and legacy coordinates:
`cos(SZA)/cos(DZA)**viewing_exponent`. Exponent zero equals `cos_sza`; exponent
one equals the legacy coordinate. It uses the same one-nightside plus
14-dayside `dgsigma` structure. Five focused tests passed, and the p=0.5 run
wrote all 129 comparison figures to
`figures/IRLS_iterations_powered_legacy_p0p5/`.

An unverified `geometry_model="powered_asymmetric_legacy"` is also available.
It first applies the powered DZA correction and then compresses the negative
coordinate using the existing `nightside_width` transformation. The code and
comparison-script label exist, but it has deliberately not been run or tested.

Two unverified spherical legacy comparisons are available. Both multiply
`cos(SZA)` by a finite spherical-shell viewing factor and match the legacy
dayside at nadir without `viewing_exponent` or `twilight_width`.
`spherical_legacy` leaves the negative coordinate unchanged, while
`spherical_asymmetric_legacy` applies the existing `nightside_width`
compression. Code, documentation, and unexecuted tests are present; nothing
has been run.

An optional spatial monotonicity constraint is now implemented for the BS
mean background. Calling `backgroundmodel_BS(..., monotonic=True)` prevents
the fitted background from decreasing with the dayglow geometry coordinate.
It does this through non-negative increments between adjacent spatial
B-spline control coefficients, independently for each temporal basis
component; temporal knots, order, and freedom are unchanged. The default is
still the legacy unconstrained fit. The first trust-region bounded solver
exhausted its iteration limit on a real sequence and was replaced by SciPy's
active-set BVLS solver. The low-cost default,
`monotonic_solver="normal"`, applies BVLS after reducing all pixels to the
small normal system. `monotonic_solver="full"` retains BVLS on the complete
weighted design as the more numerically stable comparison. Conditional
uncertainty is rejected for this option because the current covariance
assumes an unconstrained linear estimator. Focused verification passed (30
tests, plus a 3-test targeted
rerun), and the synthetic fitted background was non-decreasing in every
frame. All six comparison configurations completed on a 20-frame WIC
sequence; their outer IRLS fits converged in two to four iterations. On one
representative variance-weighted configuration, normal/full took 1.630/7.958
s and differed in `dgmodel` by at most (1.52\times10^{-8}).

The experimental background-fit uncertainty architecture through Phase 5 is
implemented and has representative complete-orbit validation. Default BS/SH
outputs remain unchanged. Opt-in diagnostics, fixed-weight conditional
coefficient covariance, corrected-image covariance, coherent coefficient
sampling, repeated BS and dependent BS-to-SH fits, compact temporal
covariance, and NetCDF metadata are available. Eighteen tests passed before
the final temporal-review corrections. Strengthened temporal tests were added
afterward but not rerun at the user's request.

Complete-orbit checks on bright orbit 214 and moderate orbit 1926 show that
fixed-weight covariance usually captures marginal fitted-background spread
and short-lag temporal correlation under a stated perturbation model. It does
not capture the repeated robust estimator's possible mean shift. The full WIC
chain is SH-dominated and has spatially varying BS--SH dependence, confirming
that its two component variances must not be added independently.

## Next action

Compare the legacy and effective-illumination coordinates on representative
complete orbits using identical non-Gaussian spatial weighting. Inspect the
background-versus-coordinate curve, held-out residual dependence on both SZA
and DZA, and the terminator artifact before tuning emission height, viewing
exponent, or twilight width.

For the eventual icBuilder handoff, use the final corrected field (`dgimg` for
SI or `shimg` for the full WIC correction) together with `img_variance` from
the measurement-variance loader. Do not request `uncertainty=True` or
propagate the explicit fitted-background covariance. This deliberately treats
the fitted background as fixed; background-model inadequacy remains a
separate systematic and the repeated-fit functions remain validation tools.

Compare `monotonic=False` and `monotonic=True` on representative complete
orbits, especially the background-versus-geometry curve and the terminator
artifact. This directly tests whether unconstrained spatial spline
oscillation caused the observed rise-and-fall behavior; it does not alter the
temporal model.

Compare the legacy and new variance-weighted BS fits on representative full
orbits, including damping sensitivity. The call is
`backgroundmodel_BS(..., input_variance="img_variance",
variance_weighted=True, tukey_scale="empirical", uncertainty=True)`. The
weighted solve uses \(1/\sqrt{v_i}\) row precision while the recommended WIC
test uses the empirical `dgsigma` scale for Tukey rejection. The existing
measurement-standardized option remains available with
`tukey_scale="measurement"`; `tukey_weights=False` disables Tukey weighting.
Legacy behavior remains the default. Use a positive Poisson offset when
loading measurement variance; zero variance at a fitted pixel is rejected.
Do not begin downstream integration until the WIC/SI13 variance assumptions
and the damping scale have been assessed.

`poisson_offset=0.5` now also applies to WIC equivalent output quanta. In the
reference WIC model it adds a 56.7 AD-squared variance floor (minimum sigma
about 7.5 AD). The default remains zero for backward compatibility.

The BS model also exposes `sza_max`; for example, `sza_max=90` makes every
nightside value use SZA 90 degrees in the spline geometry and empirical
`dgsigma` calculation. It does not alter `imgs["sza"]`. No clipping is the
default.

Audit of constant WIC `dgsigma` found a mismatch between the publication and
released `_noisemodel()`: a bin-level RMSE is duplicated onto all pixels and
then fitted against pixel-level background values. On all three example WIC
sequences this forced the fitted variance slope effectively to zero. The new
`noise_model="binned"` option uses paired weighted means of residual squared
and model squared, restoring nonconstant sigma in both 20-frame WIC examples.
`noise_model="legacy"` remains the default. Do not change the single nightside
bin or nonnegative coefficients; both are stated in the paper.

The next background-method question is the package's so-called normal spatial
weight. It is an inverse assumed-normal-density multiplier rather than the
paper's inverse empirical bin-population weight. It can strongly upweight
large `abs(cos(SZA)/cos(DZA))`.

The spatial comparison is now selectable through `spatial_weighting`:
`"gaussian"` retains the package default, `"bin_population"` gives each
populated geometry bin equal total least-squares weight, and `"none"` uses
unit multipliers. Bin-population row weights use inverse square-root counts
and are normalized to preserve their overall squared scale; using literal
inverse counts caused the fixed damping to collapse a representative WIC
background almost to zero. The old `gaussian_weights` input remains for
compatibility but cannot be combined with `spatial_weighting`.

The existing `figures/IRLS_iterations/` output also used the misspelled
`gaussian_weight=False`; the argument was silently ignored by the `**kwargs`
interface. All those freeze-iteration figures therefore retain the default
inverse-Gaussian weighting. They nevertheless show a repeatable
later-iteration Tukey-weight structure near
`cos(SZA)/cos(DZA)=0`, sometimes accompanied by a corrected-image artifact.
Use `spatial_weighting="none"` when regenerating them.

The zero location has a mixed physical and algorithmic interpretation. With
the example DZA limit below 90 degrees, fraction zero is exactly the solar
terminator (SZA 90 degrees), where a one-coordinate dayglow model is most
likely to leave coherent residual structure. `dgsigma` follows the dayglow
gradient by construction because it is evaluated as
`sqrt(a * dgmodel**2 + b)`. The noise-model statistics also treat zero as a
hard boundary: all negative fractions form one nightside bin while positive
fractions use 0.25-wide bins. That asymmetry can influence the fitted scale,
but cannot by itself make a discontinuous `dgsigma` field. The sharp Tukey
feature is most plausibly produced when the smooth residual-to-scale ratio
crosses the Tukey rejection threshold and later scale recalculation feeds that
rejection back into the background fit. A simple cubic-spline knot at zero
does not permit a jump, and the inverse-Gaussian multiplier is smooth there.

Two controls now support that comparison: `gaussian_weights=False` removes the
inverse-Gaussian row multiplier, and `dgsigma_freeze_iteration` selects the last
zero-based iteration that recalculates `dgsigma`. Zero is the conventional
initial-fit freeze; one is the unchanged fuvpy default. Seven focused noise,
geometry, default-reference, and control tests passed.

The fit can now retain candidate scale evolution independently of the chosen
freeze. With diagnostics enabled, `dgsigma_diagnostic_iterations=5` stores
`dgsigma_iteration` and `dgweight_candidate_iteration` for the first five IRLS
iterations (or fewer if the fit converges first). Later candidate fields are
not used by the operational fit. This should now be used to determine whether
the near-constant `dgsigma` seen in several comparison cases is present from
the initial fit or emerges during robust reweighting. The focused suite now
passes 26 tests.

For calibration-flag-1 SI12 files without FUVVIEW3 airglow subtraction, the
saved image is (y_i=a_iN_i). Here (N_i) is an integer photon-event count and
(a_i) is the date- and pixel-dependent product of mission-gain and flat-field
corrections. The physical variance is (a_i^2\lambda_i), with plug-in estimate
(a_i y_i). This relation recovered integer detector counts from every
positive usable pixel in the 40-frame 2000-day-137 example sequence, with
maximum numerical error (4.7\times10^{-6}) count.

The FUVVIEW3 mission and flat-field tables are now packaged unchanged under
`fuvpy/data/calibration/`. `read_idl(..., measurement_variance=True)` uses them
to add the combined multiplier and detector-counting variance before optional
WIC reflattening. The loader rejects unsupported SI calibration flags,
nonzero FUVVIEW3 airglow subtraction, and missing mission dates rather than
silently applying the simple equations. Ordinary calls remain unchanged.

The background-uncertainty go/no-go test remains stopped at Step 3 for the
three-sensor empirical model. Its executable baseline is recorded in
[[Background Uncertainty Frozen Mean Pipeline]], fitted models in [[Effective
Observation Variance - Step 2 Results]], and held-out validation in [[Effective
Observation Variance - Step 3 Results]].

The freeze exposed two downstream constraints. The accepted detector-first
architecture is not yet an executable pipeline, so the first screen uses the
current bin-first path. In addition, the current 46-by-46 WIC grid does not
match the bundled 36-by-36 Zhang--Paxton lookup. Resolve that mismatch before
Step 6; it does not block Step 2.

### Step 2 model and Step 3 result

Ten well-separated common orbits produced provisional effective variance laws
in corrected-count space:

- WIC: \(\sigma^2=19.0424B\), weighted \(R^2=0.944\);
- SI12: \(\sigma^2=0.3317+1.3407B\), weighted \(R^2=0.984\); and
- SI13: \(\sigma^2=2.5569B\), weighted \(R^2=0.850\).

The fits are candidates rather than accepted models. WIC's dynamic range
inflates the apparent fit quality, SI13 shows high-background curvature, and
WIC's orbit-level residual center and spread vary substantially. Residuals are
also spatially correlated: median WIC correlation is 0.75 at one detector
pixel and 0.46 at five pixels. SI correlations are weaker but persist to five
pixels. Independent-pixel perturbations are therefore not defensible for WIC.

Held-out orbits 480, 800, 1200, and 1450 gave pooled standardized spreads of
0.611 for WIC, 1.010 for SI12, and 0.645 for SI13. SI12 passed every screening
check. WIC and SI13 failed orbit-transferability and background/SZA trend
checks; SI13 additionally failed the DZA trend check. Background-only variance
is therefore not transferable for the required three-sensor system.

The next user decision is whether to record this route as closed or start a
new model-development cycle with additional predictors such as SZA and mission
time. The latter is scientifically possible but materially broader than the
bounded go/no-go test.

## Important evidence

- Count-scale adjacent-frame correlation, fixed/repeated:
  SI13 0.979/0.978; WIC BS 0.980/0.984.
- Count-scale repeated/fixed BS sigma: SI13 1.149 median; WIC 1.142 median.
- Corrected-image fixed one-sigma coverage: 0.70 median in both bright cases.
- Count-scale repeated-mean shift: -0.095 SI13 counts; +0.730 WIC counts.
- Full WIC component sigma: BS 0.238, SH 4.283, combined 4.301 counts.
- Full WIC pixelwise BS/SH correlation: median 0.025,
  16th--84th percentiles -0.253--0.337.
- Full WIC combined frame-mean adjacent correlation: 0.995.

Ten to twenty repeated fits resolve marginal and short-lag behavior but not a
publication-quality long-lag covariance matrix.

## Risks and boundaries

- SI12 has a viable detector-count model and WIC has a first-order
  equivalent-quanta candidate; SI13 remains unresolved.
- The WIC candidate assumes `IMAGE` follows the flag-1 corrected-AD path. The
  recovered source lacks the producer-side flag-3 branch, and the published
  WIC noise is specified per 2-by-2-pixel resolution cell rather than as a
  full detector-pixel covariance.
- Conditional covariance freezes basis, mask, damping, and converged weights.
- Background-model discrepancy is excluded by the stated assumption that the
  selected background model is correct.
- Calibration, flat-field, atmospheric-response, and regularization bias are
  outside the reported uncertainty.
- Public field names remain experimental.
- `polplot` is imported by the source but is not declared in `pyproject.toml`.

## Verification entry points

- `fuvpy/src/background.py`
- `fuvpy/src/background_fit_outputs.py`
- `fuvpy/src/background_uncertainty.py`
- `fuvpy/src/measurement_uncertainty.py`
- `vault/02_Algorithm/Background Uncertainty Go-No-Go Test Plan.md`
- `vault/02_Algorithm/IMAGE FUV Measurement Uncertainty - SI and WIC.md`
- `tests/test_background_uncertainty.py`
- `examples/scripts/background_uncertainty_validation.py`
- `vault/02_Algorithm/Background Fit Uncertainty Architecture.md`

Exploratory complete-orbit figures are under `/tmp/fuvpy_validation/` and are
not tracked. No commit, push, icBuilder change, or production-data rewrite was
performed.

The sensor-variance calculators, calibration-table reconstruction, opt-in
loader path, and focused formula tests were added on 2026-09-03. Syntax and
`git diff --check` passed. One-frame SI12 and WIC smoke checks confirmed the
expected shapes and formulas; WIC was checked with reflattening disabled so
the returned image remained the image used by the variance calculation. The
pytest suite was not executed.

All 43 example frames for each sensor were then processed in three contiguous
time segments. The ordinary and uncertainty-enabled BS products agree in the
inspected frame comparisons, including model, corrected image, and weights.
The run generated 129 figures under
`figures/measurement_variance_bs_comparison/`. WIC was loaded with
`reflat=False`. During the run, post-boom zero multipliers exposed an IDL
rotation mismatch; matching IDL's nearest-edge extrapolation resolved it.

The BS fit now also supports explicit opt-in variance-weighted IRLS. This is a
mean-fit change, unlike the earlier post-fit covariance calculation. Five
focused tests passed, covering both Tukey scales, the disabled-Tukey path,
conditional covariance, and rejection of zero fit variance. Python compilation and
`git diff --check` passed. The comparison script was updated, but the existing
129 figures were not regenerated and still show the earlier post-fit-only
comparison.

## Portfolio impact

- Central update needed: Yes
- Changes: The physical SI and candidate WIC detector-variance calculations
  are now available as an opt-in loader product using packaged FUVVIEW3
  calibration tables. SI13 and WIC still require scientific validation, so
  nothing should yet be propagated into icBuilder.
- Sync summary: fuvpy can now produce `img_variance` at load time, but the
  required three-sensor model is not yet accepted.
