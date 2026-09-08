# Background Fit Uncertainty Architecture

Last reviewed: 2026-09-07
Status: experimental implementation and complete-orbit validation through
Phase 5; detector observation variance remains unresolved

## Implementation checkpoint

The fixed-weight BS and conditional SH covariance, corrected-image diagonal,
coefficient-space sampling, full repeated BS and BS-to-SH fits, diagnostics,
compact temporal covariance, and NetCDF fields are implemented on
`background_uncertainty`. Synthetic, serialization, and complete-orbit tests
pass. The implementation remains experimental because no calibrated detector
observation-variance model has yet been selected.

## 1. Objective

Estimate the uncertainty introduced when fuvpy fits and subtracts image
backgrounds. The work must preserve the existing fitted background and
corrected image when uncertainty calculation is disabled.

The first useful result is not a claim of total image uncertainty. It is a
tested estimate of **conditional fitted-background uncertainty**: how uncertain
the fitted background is under a stated observation-noise model, conditional
on the selected basis, regularization, fit mask, and converged robust weights.

This estimate can later become one component of detector-level uncertainty in
icBuilder. It must not be presented as including calibration, flat-field,
atmospheric-response, or wrong-background-model errors.

## 2. Scope

Initial scope:

- the B-spline dayglow fit in `backgroundmodel_BS()`;
- the spherical-harmonic residual fit in `backgroundmodel_SH()`;
- the dependence of the WIC SH fit on the preceding BS fit;
- analytic conditional covariance and repeated-fit validation; and
- compact outputs that retain spatial and temporal correlation.

Outside the initial scope:

- replacing the mean background model (an opt-in monotonicity experiment is
  documented below);
- selecting new knots, damping, or Tukey constants;
- camera calibration or reflattening uncertainty;
- atmospheric model discrepancy;
- the alternative `backgroundmodel_SBS()` implementation;
- changes to icBuilder products; and
- a dense pixel-by-pixel covariance matrix.

## 3. Uncertainty quantities

The implementation must keep the following quantities distinct.

### 3.1 Observation or count uncertainty

This describes uncertainty in the input detector values before background
fitting. It may eventually include count statistics and other quantified
detector effects. It should enter an analytic propagation as an explicit input
variance, not be inferred silently from the fit weights.

### 3.2 Residual spread

The existing `dgsigma` models the empirical spread of residuals as a function
of fitted signal. It is used to normalize the SH calculation. The residuals
can contain aurora, detector noise, imperfect background structure, and
outliers. Therefore:

> `dgsigma` is not the standard deviation of `dgmodel` and is not a total
> uncertainty on `dgimg`.

It can be tested as a provisional observation-scale model, but that use must
be named explicitly.

### 3.3 Conditional fitted-background uncertainty

This is uncertainty in fitted coefficients and predicted background after
freezing:

- the basis and knots;
- the fit mask;
- the damping matrix;
- the converged Tukey weights;
- the spatial sampling weights; and
- any supplied observation variances.

This is the first analytic target.

### 3.4 Background-model discrepancy

This is error caused by the chosen BS or SH model being unable to represent
the true background. It also includes omitted physical dependencies. A
coefficient covariance around the selected fit does not measure this error.
Repeated fits to synthetic noise do not measure it either unless the model
itself is varied.

## 4. Existing calculations

### 4.1 B-spline fit

For the selected pixels, let:

- `d` be the flattened detector values;
- `G` be the joint spatial and temporal B-spline design matrix;
- `w` be the converged Tukey weights;
- `ws` be the spatial weights; and
- `R` be the diagonal damping matrix.

The current code multiplies each row by

```text
q = w * ws
```

before constructing the normal equations. It therefore solves

```text
H m = G.T Q^2 d

H = G.T Q^2 G + R
Q = diag(q)
```

or, equivalently,

```text
minimize || Q (d - G m) ||^2 + m.T R m.
```

The square on the combined weights is important. An uncertainty refactor must
first reproduce this exact calculation, whether or not a different weighting
convention might be preferable scientifically.

The code currently retains the fitted image, corrected image, final Tukey
weights, residual-spread model, and two norms. It does not retain the fitted
coefficient vector, normal matrix, coefficient ordering, convergence state, or
matrix-conditioning diagnostics.

### 4.2 Noise model

`_noisemodel()` bins residual RMS values and fits

```text
residual_variance = a * fitted_signal^2 + b.
```

SciPy returns a two-by-two covariance for `a` and `b`, but it is discarded.
That small covariance describes the fitted residual-spread relation. It is not
the covariance of the BS background coefficients. Moreover, `_noisemodel()`
is evaluated only during IRLS iterations zero and one. The resulting
`dgsigma` is then frozen while the robust fit continues, so it describes an
early-fit residual spread rather than the residuals of the converged fit.

The discarded two-parameter `curve_fit` covariance is itself only a formal
fit covariance. The same binned RMSE is repeated for every pixel in a bin, so
it must not be treated as a calibrated uncertainty on `a` and `b` without a
separate statistical justification.

The experimental implementation now separates the operational scale from its
diagnostic evolution. `dgsigma_freeze_iteration` controls the last candidate
used by IRLS. When diagnostics are enabled,
`dgsigma_diagnostic_iterations` independently controls how many early
candidate scale fields and corresponding empirical Tukey-weight fields are
retained. Candidates calculated after the freeze are not fed back into the
fit. This makes it possible to diagnose whether a constant fitted scale is a
property of the binned variance law or a later consequence of reweighting.

The BS spatial sampling multiplier is independently selectable through
`spatial_weighting`. `"gaussian"` retains the package behavior,
`"bin_population"` gives every populated geometry bin equal total
least-squares weight, and `"none"` uses unit multipliers. Because row
multipliers are squared by the normal equations, bin-population weighting uses
inverse square-root counts. It is normalized to retain unit mean squared
weight so that selecting it does not inadvertently change the scale of the
fixed damping.

### 4.3 Spherical-harmonic fit

The SH stage fits `dgimg / dgsigma`. Its design matrix is also divided by
`dgsigma`. It starts with the BS Tukey weights and adds equal-area spatial
weights before another robust iterative solve.

Any conditional SH covariance must state its units. The coefficient fit is in
the normalized system, whereas `shmodel` is converted back to detector units
by multiplying its prediction by `dgsigma`. The prediction covariance must be
transformed back in the same way. Such a calculation also freezes `dgimg`,
`dgsigma`, and the inherited BS weights.

The final WIC correction is

```text
shimg = img - dgmodel - shmodel.
```

Consequently, BS and SH errors are dependent. The SH fit changes when the BS
fit, `dgsigma`, or BS weights change. Separate BS and SH variances must not be
added as if the stages were independent.

### 4.4 Optional spatial monotonicity experiment

The unconstrained spatial B-spline can rise and fall with the geometry
coordinate even without IRLS. Adding more knots increases this freedom and
does not enforce the physically expected behavior. An experimental
`monotonic=True` option therefore constrains adjacent spatial control
coefficients to be non-decreasing.

The coefficient array is spatial by temporal. The constraint is applied
independently down the spatial direction for every temporal basis component.
It does not change the temporal knots, order, or coefficient freedom. Because
B-spline basis functions are non-negative, this construction makes the
background non-decreasing with the spatial geometry coordinate at every
modeled time.

The constraint is opt-in and the legacy fit remains the default. The existing
conditional covariance is not used with this option: active inequality
constraints make the estimator piecewise linear rather than the
unconstrained linear smoother assumed by that covariance calculation.

Two bounded solvers are retained. `monotonic_solver="full"` applies BVLS to
the complete weighted pixel design and is the more numerically stable
reference. `monotonic_solver="normal"` first forms the same small normal
matrix and right-hand side used by the original fit, transforms that system
to the monotonic-increment parameters, and applies BVLS only to the resulting
coefficient-sized problem. The latter is the low-cost default, with the usual
normal-equation tradeoff of squaring the design matrix's condition number.

### 4.5 Optional effective-illumination coordinate

The legacy coordinate `cos(SZA)/cos(DZA)` maps every terminator pixel to zero
and cannot distinguish viewing-path differences there. The experimental
`geometry_model="effective_illumination"` instead uses

```text
x = log(smooth_solar(SZA)) + alpha * log(van_Rhijn(DZA)).
```

The stored coordinate is translated by its value at SZA 90 degrees and DZA 0
degrees. The nadir terminator is therefore zero. Default knots receive the
same translation, so this choice of origin does not change the spline fit.

The smooth solar term approximates positive `cos(SZA)` without a hard break
at 90 degrees. The van Rhijn term represents the path length through a thin
spherical emitting shell. The default parameters are a 130 km shell height,
`alpha=1`, and a 5 degree twilight width. The logarithm compresses the
day-to-night range without changing the coordinate ordering.

The coordinate remains one-dimensional and can use the same monotonic BS
solver. Its default knots are fixed by interpretable SZA locations with extra
resolution around the terminator, plus a knot at unit effective illumination
to split the otherwise broad bright-side interval.
The legacy inverse-Gaussian sampling weight is not transferred because its
center and scale were defined for the old coordinate; neutral weighting is the
new default and bin-population weighting remains available.

The alternative
`geometry_model="bounded_effective_illumination"` uses the same effective
illumination, `I`, but transforms it as

```text
x = log(1 + I / I_terminator) - log(2).
```

Here `I_terminator` is the smooth-solar illumination at SZA 90 degrees and
DZA 0 degrees. This retains zero at the nadir terminator while bounding the
deep nightside at `-log(2)`. It avoids the raw logarithm's strong expansion of
very small nightside values. The raw-log and bounded versions are separate
options; neither replaces the legacy default.

The additional `geometry_model="linear_effective_illumination"` uses

```text
x = (I - I_terminator) / (I_subsolar - I_terminator).
```

It preserves the effective-illumination geometry without logarithmic
compression. At nadir the terminator is zero and the subsolar point is one.
Default knots include nadir SZA values of 75, 60, 45, and 30 degrees to
resolve the dayside gradient, followed by oblique dayside anchors at DZA 60,
75, and 90 degrees.

The separate `geometry_model="asymmetric_legacy"` starts from

```text
u = cos(SZA) / cos(DZA)
```

and uses

```text
x = u                                      for u >= 0
x = nightside_width * tanh(u / width)      for u < 0.
```

The dayside is therefore identical to the legacy coordinate. The nightside
retains the same value and first derivative at the terminator but rapidly
saturates at `-nightside_width`; the default width is 0.2. The default
negative knots scale with this width, while the positive legacy knots remain
0.25, 1.5, and 3.5.

## 5. Conditional analytic covariance

Assume for the moment that the converged weights and all fit choices are
fixed. Let the input measurement covariance be `Sigma_d`.

When `H` is nonsingular, the fitted coefficient estimator can be written as

```text
m_hat = H^-1 G.T Q^2 d.
```

and its frequentist covariance under the stated input covariance is

```text
C_m = H^-1 G.T Q^2 Sigma_d Q^2 G H^-T.
```

The current code permits zero damping and solves the normal system by least
squares. `H` is therefore not guaranteed to be nonsingular. The implementation
should reproduce the legacy generalized solution and define

```text
K = solve_legacy(H, G.T Q^2)
C_m = K Sigma_d K.T.
```

The inverse notation above is only shorthand for the full-rank case. The code
must use linear solves or a factorization, not an explicit inverse, and report
both the rank of the weighted design `Q G` and the conditioning of `H`.

For prediction design matrix `G_pred`, the fitted-background covariance is

```text
C_background = G_pred C_m G_pred.T.
```

Only its diagonal is required for a per-pixel standard-deviation field:

```text
var_background[i] = G_pred[i] C_m G_pred[i].T.
```

The diagonal can be calculated row by row without constructing the dense
pixel covariance.

This is the sampling covariance of the fitted background estimator. It does
not include regularization bias. If damping materially shrinks the fit, that
bias must be assessed separately rather than folded into `C_m`.

### Corrected-image covariance

The fitted background and input detector image are derived from the same
data. Their errors are therefore correlated. With a fixed linear smoother
`S`,

```text
background_hat = S d
corrected_image = (I - S) d
```

and

```text
Cov(corrected_image) = (I - S) Sigma_d (I - S).T.
```

It is generally wrong to add the original detector variance and
`C_background` as independent terms. Repeated fits preserve this dependence
when both the perturbed input and its fitted background are carried through.
Drawing coefficient perturbations around a fixed, unperturbed detector image
does not.

### Why `H^-1` alone is not the answer

`H^-1` could be interpreted as a Bayesian posterior covariance only if `Q^2`
is a calibrated observation precision and `R` is a calibrated Gaussian-prior
precision. The current Tukey and spatial weights were not defined that way,
and default damping may be zero. The initial implementation must therefore
call the result a conditional covariance, not a posterior covariance.

Likewise, Gaussian coefficient draws using `C_m` are a practical approximation
to the sampling distribution of the fitted estimator. They are not posterior
draws unless an explicit likelihood and prior have first been defined.

## 6. Proposed software shape

Keep the change small and readable.

### 6.1 Freeze reference behavior first

Before moving any solve code, save reference outputs for a small deterministic
dataset. Tests should cover `dgmodel`, `dgimg`, `dgweight`, `dgsigma`, the
norms, and the equivalent SH fields.

### 6.2 Expose fit diagnostics internally

Extract one short internal weighted-solve helper only after the reference test
exists. It should return ordinary arrays in a dictionary, for example:

```text
coefficients
normal_matrix
combined_row_weights
fit_indices
rank
condition_number
```

The main BS and SH workflows should remain visible from top to bottom. The
helper should not become a generic fitting framework.

The diagnostic result must also state the coefficient order. In the present
BS code the temporal coefficient index varies fastest within each spatial
basis coefficient.

### 6.3 Add covariance as an opt-in calculation

The public default remains unchanged and performs no uncertainty calculation.
The experimental interface is:

```python
backgroundmodel_BS(imgs, uncertainty=True, ...)
```

The same opt-in convention is available for SH. The names remain experimental
until complete-orbit evidence exists.

### 6.4 Use explicit output names

The implemented per-pixel names are:

```text
dgmodel_fit_sigma
shmodel_fit_sigma
```

These names deliberately say that the values describe fitted-model
uncertainty. Generic names such as `sigma` or `image_uncertainty` would invite
misinterpretation.

## 7. Observation variance supplied to the calculation

The covariance routine should accept a one-dimensional variance array for the
selected detector samples. It should not infer that variance from Tukey or
spatial weights.

Three inputs should be compared during development:

1. known variances in synthetic tests;
2. a provisional variance based on `dgsigma^2`; and
3. a later detector/count variance supplied by the data pipeline.

The second option is useful diagnostically but must be labelled
`early_fit_residual_scale` rather than detector noise. The comparison will
show how much the uncertainty estimate depends on this choice.

## 8. Data-dependent robust weights

The analytic equation freezes the final Tukey weights. In reality, perturbing
the image can change which pixels are downweighted. Repeated complete fits are
therefore needed for validation. A robust sandwich estimate is an optional
intermediate comparison, not a required product.

### 8.1 Robust sandwich estimate

A sandwich or M-estimator covariance may account approximately for the robust
score and the observed residual distribution. Implement it only if the
fixed-weight result disagrees enough with repeated fits to justify the added
complexity. The present code's unusual row-weight convention and early-frozen
scale model would have to be represented exactly.

### 8.2 Complete repeated fits

Repeated fits are the Monte-Carlo reference under the selected perturbation
model. Each realization should perturb the input sequence, then rerun the
complete IRLS calculation so the robust weights and `dgsigma` can change. The
ordinary additive perturbations proposed here do not change the geometric and
finite-data fit mask; missingness must be perturbed deliberately if mask
uncertainty is to be studied.

Candidate perturbations are:

- parametric draws from an explicit heteroscedastic observation model;
- wild-bootstrap residuals; and
- spatial or temporal block residuals if residual correlation is important.

Naively resampling all residual pixels independently is not sufficient: the
residual field contains auroral structures and correlations that the fit is
designed to reject.

The repeated-fit ensemble must use a complete orbit or complete fitted image
sequence. The temporal spline problem changes when only a few extracted
frames are supplied.

## 9. Sequential WIC uncertainty

icBuilder currently consumes the BS-corrected `dgimg` field for SI12 and SI13,
so the BS calculation is the natural first target. fuvpy itself permits either
background routine to be applied to any supported sensor.

WIC normally uses the output of both BS and SH stages. Development should
proceed in this order:

1. calculate BS uncertainty with SH disabled;
2. calculate SH uncertainty conditional on the fitted BS result;
3. rerun the entire BS-to-SH chain for every repeated-fit realization; and
4. compare the combined WIC spread with the two conditional analytic pieces.

The complete repeated-fit WIC spread is the reference for the combined
background-fit uncertainty. Adding `dgmodel_fit_sigma^2` and
`shmodel_fit_sigma^2` is explicitly rejected because it omits their covariance.

A joint fixed-weight linear propagation through both stages may later be
possible, but it is not required for the first validated result.

## 10. Preserve correlation compactly

A pixelwise standard deviation is useful for plots and a heteroscedastic VAE
likelihood, but it does not preserve the correlation introduced by the
low-dimensional background fit.

The preferred representation is coefficient space:

- fitted coefficient vector;
- coefficient covariance or a Cholesky/eigen factor;
- coefficient ordering;
- exact spatial and temporal knots and orders;
- temporal origin and frame times;
- prediction and mask semantics;
- for SH, the harmonic keys and frozen `dgsigma`; and
- enough metadata to reconstruct the prediction design.

This is small compared with a dense pixel covariance. A downstream sampler can
draw coefficient perturbations and map them through the design matrix to
generate spatially and temporally coherent background realizations.

These arrays are stored inside the xarray image Dataset. Synthetic and bundled
frame tests confirm a NetCDF round trip and refuse sampling when current code
would reconstruct a basis inconsistent with the saved knots or harmonic keys.

### Complete-orbit validation outcome

Validation used complete detector-grid sequences rather than extracted
frames. Orbit 214 provided a bright 57-frame WIC case and a 56-frame SI13
case; orbit 1926 provided moderate 64-frame WIC and 65-frame SI13 cases.

Two deliberately different perturbation scales were tested:

- the provisional early-fit residual scale, ``dgsigma**2``; and
- ``max(img, 1)``, labelled only as a Poisson/count-scale diagnostic.

Across these cases, repeated/fixed marginal background-sigma ratios were
typically about 1.0--1.15. Corrected-image one-sigma coverage had a median of
about 0.70, close to the Gaussian reference value. Short-lag temporal
correlation was also reproduced: for the count-scale bright case, analytic
and repeated adjacent-frame correlations were 0.979/0.978 for SI13 and
0.980/0.984 for WIC.

The repeated fitted background was not always centered on the unperturbed
robust fit. Under the count-scale diagnostic, median absolute shifts were
about -0.095 SI13 counts and +0.730 WIC counts. These are small in native
units but can be several conditional fit sigmas because that conditional
sigma is very small. A covariance cannot represent this nonlinear estimator
mean shift.

The residual-scale model is not a defensible detector-noise default. In the
tested WIC sequences, ``dgsigma`` was roughly 3--5.3 times ``sqrt(img)`` at
the median, so it also contains residual auroral/model structure. It produced
much larger nonlinear mean shifts.

For the full count-scale WIC BS-to-SH chain, the median repeated component
sigmas were 0.238 and 4.283 native units for BS and SH, respectively; the
combined value was 4.301. Pixelwise BS/SH correlation changed sign across the
detector, with median 0.025 and 16th--84th percentiles -0.253--0.337. The
combined frame-mean background had adjacent-frame correlation 0.995. This
confirms that the dependent repeated chain, not an independent variance sum,
is the relevant WIC reference.

These results support the fixed-weight covariance as a compact description of
conditional marginal spread and short-lag temporal correlation for a stated
variance model. They do not validate the observation variance, remove a
nonlinear fit mean shift, quantify model discrepancy, or define total image
uncertainty. Long-lag repeated correlation from only 10--20 realizations is
too noisy for a scientific covariance estimate. The implemented temporal
summary is an unweighted mean over valid detector pixels; it is a diagnostic,
not an area-weighted polar mean. Its lag plot uses frame-index lag, which is
equivalent to elapsed time only for a sufficiently regular cadence.

## 11. Implementation phases

### Phase 0 — Reference and terminology

- Add deterministic BS and SH reference tests.
- Record exact output fields, shapes, coefficient order, and weight semantics.
- Add convergence, `rank(Q G)`, and conditioning diagnostics without changing
  the fitted result.

**Exit criterion:** uncertainty-disabled output matches the current code to a
strict numerical tolerance.

### Phase 1 — Conditional BS covariance

- Retain the final BS coefficient vector and normal matrix internally.
- Implement covariance propagation from an explicit input variance.
- Calculate the diagonal fitted-background prediction variance efficiently.
- Verify it on a synthetic linear problem with known covariance.
- Verify coefficient draws reproduce the analytic prediction variance.

**Exit criterion:** the covariance is finite, symmetric, positive
semidefinite within numerical tolerance, and sampling agrees with the analytic
diagonal.

### Phase 2 — Robustness comparison

- Perturb synthetic data strongly enough to change Tukey weights.
- Compare fixed-weight and repeated-fit estimates.
- Add a robust sandwich candidate only if that comparison shows it would be
  scientifically useful.
- Document when the conditional approximation becomes overconfident.

**Exit criterion:** the selected analytic approximation has a measured and
reported range of validity.

### Phase 3 — Complete-orbit BS validation

- Select one complete SI sequence and one complete WIC sequence.
- Run repeated complete BS fits using a stated perturbation model.
- Compare spatial patterns, temporal correlation, and empirical coverage.
- Repeat with at least one quiet and one bright/contaminated case when data are
  available.

**Exit criterion:** analytic and repeated-fit differences are quantified on
representative data, not just the bundled partial example.

### Phase 4 — Sequential WIC BS and SH

- Add SH fit diagnostics and conditional covariance.
- Rerun the full BS-to-SH chain for every WIC realization.
- Estimate the combined corrected-image spread and correlation.
- Determine whether a useful compact joint representation can be constructed.

**Exit criterion:** the reported WIC value respects BS/SH dependence and is
not an independent variance sum.

### Phase 5 — Public API and serialization

- Choose stable field names and attributes with clear units and scope.
- Test xarray and NetCDF round trips for the compact representation.
- Document how to obtain a diagonal field and coherent samples.
- Keep all new work opt-in until validation is complete.

**Exit criterion:** another scientist can identify exactly what is and is not
included in every uncertainty field.

### Phase 6 — Downstream integration

- Pass the diagonal component to icBuilder only if needed for per-pixel
  likelihoods or diagnostics.
- Pass samples or a compact factor when correlation must survive Monte Carlo
  propagation.
- Keep background-fit uncertainty separate from count, calibration,
  flat-field, and model-choice uncertainty in downstream metadata.

This phase belongs to a later icBuilder task and must not be folded into the
initial fuvpy implementation.

## 12. Tests and acceptance checks

At minimum, the uncertainty implementation needs:

- legacy-equivalence tests with uncertainty disabled;
- a small ordinary least-squares case with a known covariance;
- a heteroscedastic synthetic case with supplied variances;
- symmetry and positive-semidefinite checks;
- rank-deficiency and poor-conditioning diagnostics;
- agreement between coefficient sampling and `diag(G C_m G.T)`;
- an outlier case that changes Tukey weights;
- a sequential WIC BS-to-SH repeated-fit test;
- a bundled fifteen-frame smoke test; and
- a deliberately slow complete-orbit validation outside the normal quick
  test suite.

The tests must preserve the existing masks of legacy fields exactly. New
uncertainty fields must be NaN wherever their prediction or variance is
undefined; NaNs must not silently become zero uncertainty. The uncertainty
mask need not equal the stricter set of samples used to fit the model.

## 13. Performance

The coefficient systems are small compared with the number of image pixels.
Coefficient covariance is therefore tractable. Dense pixel covariance is not.

Use:

- one factorization per converged fit where possible;
- triangular or least-squares solves rather than explicit inverses;
- row-wise prediction-variance calculations;
- coefficient-space draws; and
- a configurable, modest number of repeated fits for scientific validation.

Repeated fits are expected to be the expensive part and should never run by
default.

## 14. Decisions still required

The architecture leaves three scientific choices open until evidence exists:

1. whether the provisional analytic result should use `dgsigma^2`, a detector
   count model, or expose both explicitly;
2. whether fixed-weight or robust-sandwich covariance is adequate relative to
   repeated fits;
3. which repeated-fit perturbation best represents the background population
   without resampling aurora as noise.

These are validation questions, not API preferences. They should be resolved
from synthetic tests and complete-orbit evidence before fuvpy or icBuilder
claims a background-fit uncertainty product.
