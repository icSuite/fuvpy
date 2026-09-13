# Background Model Revision Implementation Plan

Last reviewed: 2026-09-08

## Objective

Create a clean, reviewable revision of the released fuvpy BS background model.
Only the changes supported by the uncertainty experiments are included. The
legacy spatial coordinate remains

\[
x = \frac{\cos(\mathrm{SZA})}{\cos(\mathrm{DZA})}.
\]

The experimental geometry coordinates, alternative knots, and generated
diagnostic figures are deliberately excluded. A monotonic solver is retained,
but it continues to use the legacy spatial coordinate.

## Branch strategy

Create `background_model_revision` from the released `main` branch, not from
`experimental`. Selectively bring across the project guidance, vault, and the
small pieces of tested implementation needed for each step. Do not merge or
cherry-pick the large experimental commit.

Each scientific change should be a separate commit. Its focused tests and
documentation belong in the same commit.

## Step 0: Create the clean branch

```bash
git switch main
git switch -c background_model_revision
```

Bring across `AGENTS.md`, the project vault, and this plan.

Suggested commit:

```text
Add project guidance for background-model revision
```

## Step 1: Correct the dgsigma estimation

Use one nightside bin for all \(x < 0\), followed by 14 dayside bins between
0 and 3.5 with a width of 0.25.

For each populated bin \(j\), calculate one pair:

\[
X_j = \left\langle I_{\mathrm{model}}^2 \right\rangle_j,
\qquad
Y_j = \left\langle (I-I_{\mathrm{model}})^2 \right\rangle_j.
\]

Fit the non-negative model

\[
Y_j = aX_j+b, \qquad a,b \geq 0,
\]

and evaluate the fitted spread at every pixel:

\[
d\sigma_i = \sqrt{aI_{\mathrm{model},i}^2+b}.
\]

Each populated bin must contribute once to the regression. Do not copy its
RMSE back to every pixel before fitting, because that weights the regression
by bin population.

Focused verification:

- confirm the one-nightside-plus-14-dayside bin structure;
- recover known \(a\) and \(b\) from synthetic data;
- confirm that duplicating observations within one bin does not change the
  fitted regression.

Suggested commit:

```text
Fix binned dgsigma estimation
```

## Step 2: Add bin-population spatial weighting

Divide the legacy coordinate into bins of width 0.1, independently for each
image frame. If a populated bin contains \(N_j\) pixels, use the row
multiplier

\[
w_{\mathrm{spatial},i} \propto \frac{1}{\sqrt{N_j}}.
\]

The squared least-squares objective then gives every populated bin equal total
influence. Normalize the multipliers to have mean square one over the fitted
pixels. This preserves the overall balance between the data and any damping
while changing their relative spatial influence.

Focused verification:

- all pixels in one bin receive the same multiplier;
- every populated bin has equal summed squared multiplier;
- weighting is calculated independently for every frame;
- empty bins have no effect.

Suggested commit:

```text
Add bin-population spatial weighting
```

## Step 3: Freeze dgsigma after iteration 0

Use this IRLS sequence:

1. Perform the initial background fit.
2. Estimate \(d\sigma\) from its residuals.
3. Calculate the first Tukey weights.
4. Keep that \(d\sigma\) fixed during all later IRLS iterations.

The revised calculation always freezes `dgsigma` after iteration 0. Later
experiments showed no production reason to expose alternative freeze
iterations, so the clean implementation has no freeze-iteration option.

Focused verification:

- the operational \(d\sigma\) equals the iteration-0 estimate;
- later model and Tukey-weight updates cannot change it;
- the IRLS convergence logic otherwise remains unchanged.

Suggested commit:

```text
Freeze dgsigma after the initial fit
```

## Step 4: Calculate measurement variance during IDL loading

Add `img_variance` while reading WIC, SI12, and SI13 IDL files.

For SI12 and SI13, use

\[
\operatorname{Var}(I_{\mathrm{SI}}) = aI_{\mathrm{SI}},
\]

where \(a\) is reconstructed from the date-dependent mission-gain and
flat-field tables.

For WIC, use the reference compound-Poisson approximation

\[
\operatorname{Var}(I_{\mathrm{WIC}}) = 10.65I_{\mathrm{WIC}}.
\]

Retain the optional 0.5 Poisson offset for zero and very low counts. Its use
must be explicit and recorded in the variable attributes. Store the small
FUVVIEW3 calibration tables in the package so the calculation is reproducible.

The variance is conditional on the calibration and background model being
correct. It excludes calibration uncertainty, background-model discrepancy,
and atmospheric-response uncertainty. WIC reflattening is treated as a
deterministic additive adjustment, so the detector-counting variance is
unchanged conditional on the estimated reflattening background.

Focused verification:

- reproduce known SI and WIC variance values;
- confirm the calibration multiplier selected for representative dates;
- confirm that dimensions, masks, units, and metadata follow `img`;
- reject unsupported sensors or calibration states clearly.

Suggested commit:

```text
Calculate detector measurement variance when reading IDL files
```

## Step 5: Use measurement variance in the BS fit

Apply the measurement row multiplier

\[
w_{\mathrm{measurement},i}
=
\frac{1}{\sqrt{\operatorname{Var}(I_i)}}.
\]

The fit therefore combines:

- bin-population spatial weighting;
- measurement-variance weighting;
- the existing Tukey robust weighting.

Keep their roles separate and visible in the code. Measurement variance is
used in the background-model inversion only. It must not be used when fitting
the empirical \(d\sigma\) model, and the Tukey calculation weights continue to use
the empirical \(d\sigma\).

Focused verification:

- constant variance reproduces the corresponding unweighted fit when damping
  is zero;
- increasing one observation's variance reduces its influence;
- changing measurement variance does not directly change the \(d\sigma\)
  calculation;
- invalid or non-positive variances used by the fit are handled explicitly.

If non-zero damping is used, its scale must be reassessed because absolute
measurement weights change the balance between the data and regularization.

Suggested commit:

```text
Use detector variance in the BS background fit
```

## Step 6: Add the monotonic spatial solver

Add an optional constraint that prevents the fitted background profile from
decreasing as the legacy spatial coordinate increases. Leave the temporal
B-spline unchanged.

Represent each spatial B-spline coefficient as a free first coefficient plus
non-negative increments:

\[
c_0 = q_0,
\qquad
c_k = q_0 + \sum_{j=1}^{k}q_j,
\qquad q_j \geq 0 \text{ for } j>0.
\]

Use the low-cost normal-system solver as the operational implementation. It
first reduces all image pixels to the same small normal system as the ordinary
fit and then solves the bounded coefficient problem. Retain the full-pixel
bounded solver as an explicit reference option for verification, not as the
default calculation.

The monotonic constraint must work with bin-population and measurement-
variance row weights. It does not introduce analytic background-model
covariance; the saved `img_variance` remains the conditional measurement
variance associated with the corrected image.

Focused verification:

- fitted spatial coefficients never decrease;
- the low-cost and full-pixel solvers give equivalent background maps on a
  small synthetic problem;
- disabling the constraint reproduces the ordinary solver;
- temporal coefficients and the legacy coordinate are unchanged;
- failure of the bounded solver produces a clear error rather than silently
  returning an unconverged model.

Suggested commit:

```text
Add monotonic spatial BS solver
```

## Step 7: Combined validation

Use a small number of complete WIC, SI12, and SI13 orbit sequences. Compare:

1. the released model;
2. the released model with only corrected \(d\sigma\);
3. the preceding model with bin-population weighting;
4. the preceding model with the iteration-0 freeze;
5. the preceding model with measurement-variance weighting;
6. the complete model with the monotonic constraint.

Inspect the background model, corrected image, \(d\sigma\), Tukey weights, and
convergence. Generated diagnostic figures remain local and should not be
committed in bulk.

Do not change the released package defaults until this combined comparison is
accepted. Initially expose the revised behaviour through explicit arguments,
so both the released result and the revised result remain reproducible.

## Relationship to the published implementation

This is a revised background model, not an exact reproduction of the paper's
accompanying code.

- The paper code bins residuals but copies each bin's RMSE back to its pixels
  before regression. The revision deliberately gives each bin one regression
  point.
- The paper code applies \(1/N_j\) as a row multiplier. Because the solver
  squares row multipliers, that gives total bin influence proportional to
  \(1/N_j\). The revision uses \(1/\sqrt{N_j}\), giving equal total influence
  per populated bin.
- The paper and released code recalculate \(d\sigma\) through iteration 1. The
  revision freezes it after iteration 0.
- Detector measurement-variance weighting is new.
- The monotonic spatial constraint is new.

These differences must be stated explicitly in later documentation and any
publication-method description.
