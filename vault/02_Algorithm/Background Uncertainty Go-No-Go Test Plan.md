# Background Uncertainty Go/No-Go Test Plan

Last reviewed: 2026-09-02
Status: Stopped at statistical gate; Step 3 failed on 2026-09-02

Progress records:

- [[Background Uncertainty Frozen Mean Pipeline]]
- [[Effective Observation Variance - Step 2 Results]]
- [[Effective Observation Variance - Step 3 Results]]

## Purpose

Determine whether uncertainty in the fitted fuvpy background is statistically
defensible and sufficiently important to propagate through the icBuilder data
products.

This is a bounded screening experiment. It must not alter the production
pipeline until the background uncertainty passes both:

1. a statistical-validity gate; and
2. a downstream-relevance gate.

The earliest data product in scope is the saved FUVVIEW3 IDL image. Level-1
UDF recovery and changes to FUVVIEW3 are out of scope. The quantity being
estimated is therefore an effective observation variance in post-FUVVIEW
corrected-count space, conditional on the selected background model.

## Scientific question

Does uncertainty in the fitted WIC, SI12, and SI13 backgrounds make a material
contribution to the uncertainty in IMAGE-derived precipitation and Hall and
Pedersen conductance after coregistration and final CS-grid binning?

The null result is also useful: if the contribution is negligible, the
experimental uncertainty machinery should not be added to icBuilder.

## 1. Freeze the mean-processing pipeline

**Completed 2026-09-02.** The exact executable baseline, source snapshots,
parameters, and known downstream mismatches are recorded in
[[Background Uncertainty Frozen Mean Pipeline]].

Use the same deterministic processing for the baseline and every ensemble
member:

- current fuvpy background models;
- existing WIC reflattening;
- established SI-to-WIC coregistration;
- Hardy proton correction;
- both image-ratio and Zhang-Paxton precipitation branches;
- Robinson conductance; and
- the final CS grid used by icAnalyzer.

Only the fitted background is varied. Hardy, Zhang-Paxton, Robinson,
coregistration, and grid parameters remain fixed so that the test isolates the
background contribution.

## 2. Estimate effective observation variance

**Completed 2026-09-02.** Ten well-separated complete-orbit fits produced
candidate sensor-specific marginal variance models. Strong WIC detector-space
correlation and substantial orbit-to-orbit variation were also found. The
models remain provisional until the held-out checks in Step 3. Full results
and limitations are recorded in
[[Effective Observation Variance - Step 2 Results]].

Use approximately ten complete, well-separated orbits for each of WIC, SI12,
and SI13.

For each sensor:

1. Run the ordinary background fit.
2. Calculate final residuals,

   \[
   r_i = I_{\mathrm{IDL},i} - \hat B_i.
   \]

3. Select background-dominated observations using the final Tukey weights.
4. Recenter the retained residuals robustly.
5. Estimate residual spread from the symmetric residual core, preventing the
   positive auroral tail from setting the variance.
6. Estimate spread in bins of fitted background intensity.
7. Fit the simplest adequate sensor-specific model, initially

   \[
   \sigma_{\mathrm{eff}}^2(\hat B)
   = \sigma_0^2 + \alpha \hat B.
   \]

8. Add a quadratic term only if held-out residual evidence clearly requires
   it.
9. Measure short-range residual correlation in detector coordinates.

Recovered FUVVIEW3 calibration tables and instrument literature provide
independent checks. SI has a direct photon-counting interpretation. WIC uses an
intensified CCD and requires a noise-equivalent gain, while the saved WIC
calibration flag 3 is not defined in the recovered FUVVIEW3 source. Literature
and tables therefore constrain the empirical model but do not replace it.

Do not substitute `dgsigma**2` for this estimate. `dgsigma` is an early-fit
residual scale used for robust outlier rejection and contains auroral and
unresolved residual structure.

## 3. Statistical-validity gate

**Completed 2026-09-02 — no-go.** Held-out validation on orbits 480, 800,
1200, and 1450 found pooled standardized spreads of 0.611 for WIC, 1.010 for
SI12, and 0.645 for SI13. WIC and SI13 also retained large trends with fitted
background, SZA, and orbit. SI12 passed, but the required three-sensor model
did not. Under the decision rule below, Steps 4--8 should not proceed with the
present variance laws. See [[Effective Observation Variance - Step 3 Results]].

Reserve another three or four complete orbits that were not used to estimate
the variance model.

For the held-out observations, calculate

\[
z_i = \frac{I_{\mathrm{IDL},i} - \hat B_i}
           {\sigma_{\mathrm{eff}}(\hat B_i)}.
\]

The variance model passes when:

- the robust spread of the retained standardized residuals is reasonably
  close to one;
- there is no major remaining trend with fitted background intensity;
- there is no major remaining trend with SZA, DZA, or mission time;
- the SI result is compatible with photon-counting expectations; and
- no sensor shows an unexplained systematic failure.

Failure at this stage is a no-go for the proposed uncertainty model. It does
not demonstrate that the true background uncertainty is negligible.

## 4. Select downstream evaluation orbits

Select four complete common-sensor orbits, distinct from the estimation and
validation orbits:

1. weak aurora and low counts;
2. ordinary moderate aurora;
3. bright aurora; and
4. strong dayglow or difficult viewing geometry.

Selection should be based on documented signal and geometry statistics, not
only on visually interesting examples.

## 5. Generate background ensembles

Begin with 20 ensemble members per evaluation orbit. This is adequate for a
screening estimate of marginal spread and dominant spatial structure, but not
for a publication-quality dense covariance matrix.

For every member:

1. perturb WIC, SI12, and SI13 in corrected-count space;
2. repeat each complete robust background fit;
3. repeat the dependent WIC BS-to-SH chain; and
4. retain the resulting background-corrected sensor images.

Use a fixed random seed. If detector-grid residual correlation is negligible
at scales relevant to the background basis, draw from the diagonal effective
variance model. If it is material, use small detector-space residual blocks
instead of independent perturbations.

Increase to 50 members only for borderline cases. Do not spend computation on
a larger ensemble before the relevance gate is evaluated.

## 6. Propagate through icBuilder

For every background realization, run

```text
background-corrected WIC/SI12/SI13
        -> SI-to-WIC coregistration
        -> proton correction
        -> image-ratio and Zhang-Paxton precipitation
        -> electron energy flux
        -> Hall and Pedersen conductance
        -> final CS-grid products
```

Use the deterministic baseline validity mask for continuous uncertainty
comparisons. Separately record how often ensemble members cross count or ratio
validity thresholds. This separates smooth uncertainty propagation from the
discontinuous effect of masking.

In the Zhang-Paxton branch, the background ensemble must not change electron
mean energy. Any downstream change comes through the inferred flux and
conductance. In the image-ratio branch, both the ratio-derived mean energy and
flux may change.

## 7. Quantify downstream impact

Calculate the following on both the native WIC grid and the final CS grid:

- ensemble mean minus the deterministic result;
- absolute and relative ensemble standard deviation;
- background variance relative to the existing count/binning variance;
- fraction of pixels changing validity classification;
- changes in the WIC/SI13 ratio and image-ratio electron mean energy;
- changes in electron energy flux;
- changes in Hall and Pedersen conductance;
- changes in frame-integrated Hall and Pedersen conductance; and
- the leading spatial mode of the background-induced uncertainty.

Reporting both grids is essential. A detector-grid contribution may be large
but irrelevant if it averages away during coregistration and CS-grid binning.
Conversely, a coherent fitted-background error may survive binning even when
its individual-pixel amplitude is modest.

## 8. Minimum diagnostic products

Produce four concise diagnostic products:

1. predicted and observed held-out residual spread for WIC, SI12, and SI13;
2. deterministic value, background standard deviation, and relative standard
   deviation maps for one representative event;
3. background contribution to ratio, mean energy, energy flux, Hall, and
   Pedersen results across all four evaluation orbits; and
4. the leading spatial uncertainty mode before and after CS-grid binning.

Avoid producing a large figure catalogue before the go/no-go decision.

## 9. Decision rule

### Go

Retain and integrate the background-uncertainty component if it passes the
statistical-validity gate and does at least one of the following:

- contributes at least 10% of the final variance over a meaningful fraction
  of auroral cells;
- materially changes pixel-validity classification; or
- produces coherent changes of approximately 5% or more in scientifically
  relevant conductance summaries.

### Limited go

If the effect is important only for weak-SI13 or image-ratio pixels, retain it
as ratio or quality-control uncertainty rather than as a general
Zhang-Paxton-based conductance component.

### No-go

If the statistically valid contribution remains negligible after
coregistration and CS-grid binning across all representative cases, document
the result and do not add the experimental uncertainty machinery to the
production icBuilder pipeline.

## Implementation boundary

Implement the experiment as a standalone validation workflow. Do not change
production fuvpy defaults, icBuilder file formats, or the production pipeline
until the decision gate has been evaluated.
