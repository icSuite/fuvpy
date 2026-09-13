# Effective Observation Variance — Step 2 Results

Last reviewed: 2026-09-02
Status: Step 2 complete; candidate models await held-out validation

## Question

Can the background-dominated residual core in the saved FUVVIEW3 images define
a simple sensor-specific effective observation variance for WIC, SI12, and
SI13?

The estimate is an **effective post-FUVVIEW variance in corrected-count
space**. It is not raw detector noise and, by assumption, does not include
background-model discrepancy, calibration error, or flat-field error.

## Data and calculation

Ten common, well-separated orbits were used:

`272, 385, 561, 679, 906, 1133, 1274, 1390, 1563, 1740`

Every selected sensor file contains 252--335 frames. The stored backgrounds
were fitted from those complete orbits. For this screening calculation, one
contiguous 40-frame block was sampled from each orbit. The block position was
stepped from early to late orbital phase across the ten orbits to avoid the
cost of many non-contiguous network reads.

For observations with final Tukey weight at least 0.5:

1. the residual was robustly recentered;
2. its negative half was reflected to avoid the positive auroral tail;
3. the resulting scale was estimated in 18 quantile bins of fitted background;
4. a non-negative linear model was fitted,

   \[
   \sigma_{\mathrm{eff}}^2(B)=\sigma_0^2+\alpha\max(B,0);
   \]

5. residual correlations were measured along detector rows and columns.

Three million pixels were sampled for WIC, 2.96 million for SI12, and 2.99
million for SI13.

## Candidate marginal variance models

| Sensor | Candidate variance [counts\(^2\)] | Weighted \(R^2\) |
| --- | --- | ---: |
| WIC | \(19.0424 B\) | 0.944 |
| SI12 | \(0.3317 + 1.3407 B\) | 0.984 |
| SI13 | \(2.5569 B\) | 0.850 |

These are candidate screening models, not accepted uncertainty laws. The SI12
relationship is close to linear over the sampled range. SI13 has visible
curvature at high background. WIC has a broad dynamic range, so its high
weighted \(R^2\) is dominated by bright-background bins and hides substantial
scatter at low background.

The orbit-level residual spread is also not stationary. The robust WIC scale
ranges from 18 to 197 counts and its residual center ranges from -58 to +41
counts. SI12 scales range from 2.0 to 4.0 counts. SI13 scales range from 1.3 to
5.6 counts. Held-out checks against mission time and viewing geometry are
therefore essential.

## Detector-space correlation

Median residual correlation across the ten orbits was:

| Sensor | Direction | Lag 1 | Lag 2 | Lag 3 | Lag 4 | Lag 5 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| WIC | row | 0.750 | 0.593 | 0.527 | 0.488 | 0.461 |
| WIC | column | 0.742 | 0.591 | 0.522 | 0.487 | 0.450 |
| SI12 | row | 0.134 | 0.115 | 0.105 | 0.095 | 0.087 |
| SI12 | column | 0.145 | 0.138 | 0.127 | 0.117 | 0.108 |
| SI13 | row | 0.165 | 0.157 | 0.144 | 0.133 | 0.124 |
| SI13 | column | 0.166 | 0.155 | 0.144 | 0.131 | 0.121 |

The WIC correlation is strong and persists beyond five detector pixels. SI12
and SI13 have weaker but nonzero correlation. These residual correlations may
contain unresolved spatial structure as well as observation noise; the test
cannot separate them while assuming the fitted background model is correct.

Consequently, independent pixel perturbations are not defensible for WIC. If
the variance model passes Step 3, repeated-fit ensembles should use compact
detector-space residual blocks or another explicitly correlated perturbation.

## Step 2 outcome

Step 2 has produced simple candidate marginal variance laws and demonstrated
that detector correlation cannot be ignored. It has **not** established that
the models are calibrated or transferable to other orbits.

The next gate is Step 3: test standardized held-out residuals on three or four
unused complete orbits for unit robust spread and trends with background, SZA,
DZA, and mission time. No Step 3 calculation has been started.

## Reproduction

- Historical script removed during the 2026-09-12 repository cleanup:
  `examples/scripts/effective_observation_variance.py`
- Scratch output: `/tmp/fuvpy_go_no_go/step2/`
- Main diagnostics: `effective_variance_models.png` and
  `detector_residual_correlation.png`
