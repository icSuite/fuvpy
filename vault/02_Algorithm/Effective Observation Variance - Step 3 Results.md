# Effective Observation Variance — Step 3 Results

Last reviewed: 2026-09-02
Status: Step 3 complete; statistical-validity gate failed

## Held-out design

The Step 2 variance laws were tested on four unused common-sensor orbits:

`480, 800, 1200, 1450`

These orbits were not used to estimate the models. Their stored backgrounds
were fitted from complete orbit sequences. A centered 40-frame block was read
from each file to control mounted-server I/O while retaining the full detector
geometry. Pixels required a final Tukey weight of at least 0.5.

The standardized residual was

\[
z = \frac{I_{\mathrm{IDL}}-\hat B}
         {\sigma_{\mathrm{eff}}(\hat B)}.
\]

Its center and spread were calculated with the same robust negative-half
estimator used in Step 2. The practical screening bands were fixed before the
result was inspected:

- pooled robust spread: 0.80--1.25;
- individual-orbit spread: 0.67--1.50;
- maximum/minimum binned spread ratio: at most 1.75; and
- absolute pooled center: at most 0.5.

These bands are screening criteria, not formal confidence intervals.

## Result

| Sensor | Pooled center | Pooled spread | Orbit spread range | Result |
| --- | ---: | ---: | ---: | --- |
| WIC | -0.071 | 0.611 | 0.388--0.834 | Fail |
| SI12 | -0.070 | 1.010 | 0.886--1.209 | Pass |
| SI13 | -0.161 | 0.645 | 0.385--1.082 | Fail |

All three residual centers are acceptably close to zero. The decisive failure
is the residual spread and its dependence on observing conditions.

For WIC, the binned standardized-spread ratios were 4.22 with background,
3.83 with SZA, and 1.27 with DZA. The Step 2 law generally predicts too much
spread, particularly for low-background and high-SZA observations.

For SI13, the corresponding ratios were 26.70, 15.05, and 2.12. This is not a
single multiplicative calibration error: the zero-intercept Step 2 law
underpredicts the lowest-background spread but overpredicts much of the
remaining range. The result also changes markedly between orbits.

SI12 passed every screening check. Its binned spread ratios were 1.42 with
background, 1.60 with SZA, and 1.16 with DZA. Its approximately count-dependent
variance is therefore a plausible sensor-specific empirical model, although
this alone does not validate a three-sensor background ensemble.

## Decision

The simple background-only variance model fails the statistical-validity gate
for the required WIC/SI12/SI13 system. Under the agreed go/no-go plan, the
project should not proceed to representative downstream ensembles or icBuilder
integration with these laws.

The failure does not show that fitted-background uncertainty is negligible.
It shows that fitted background intensity alone does not provide a transferable
effective variance for WIC and SI13. Adding SZA, mission-time, or orbit-specific
terms would define a new model requiring a new estimation and validation
cycle; that is a separate decision rather than Step 4 of the present test.

## Reproduction

- Historical script removed during the 2026-09-12 repository cleanup:
  `examples/scripts/validate_effective_observation_variance.py`
- Step 2 model: `/tmp/fuvpy_go_no_go/step2/variance_models.json`
- Scratch output: `/tmp/fuvpy_go_no_go/step3/`
- Diagnostics:
  - `held_out_background_spread.png`
  - `held_out_standardized_spread.png`
  - `held_out_orbit_summary.png`
