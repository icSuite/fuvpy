# FUVPY Production Reprocessing Checklist

Last reviewed: 2026-09-13
Branch: `background_model_asymmetry`

This file is the numbered source of truth for work required in fuvpy before
reprocessing the IMAGE archive. It supersedes the earlier conversational list.
icBuilder integration and restart handling are deliberately outside this list.

## 1. Remove obsolete processing paths — complete

- Remove `backgroundmodel_SH()` and `backgroundmodel_SBS()`.
- Remove code, examples, and documentation used only by those paths.
- Retain compact historical SH fields only as comparison data.
- Keep `backgroundmodel_BS()` as the sole production correction, with optional
  directional asymmetry.

## 2. Restore measurement-variance metadata — complete

`img_variance` now records:

- that it is calculated from `img` before optional WIC reflattening;
- that calibration, background-model discrepancy, and atmospheric-response
  uncertainty are excluded; and
- that WIC reflattening is treated as deterministic, so the variance is not
  changed by reflattening.

This was a metadata repair only. It did not change the calculated variance or
the fitted background.

## 3. Lock the monotonic-solver failure policy — complete

`monotonic_fallback` controls whether the ordinary solver is attempted after a
failed monotonic solve. Solver failures no longer terminate a long processing
run:

- with fallback disabled, the function returns an all-NaN background product
  with flag 3;
- with fallback enabled, a successful ordinary fit is returned with flag 2
  and `monotonic_constraint="none"`; and
- if the ordinary fallback also fails, an all-NaN product is returned with
  flag 4.

All-NaN failure products record `dgfit_converged=False`, even when the optional
diagnostics were not requested. `dgsolver_flag` stores the solver outcome and
the solver error message is retained in its metadata.

## 4. Make each BS result self-describing — complete

Save the configuration needed to understand and reproduce each background fit,
including at least:

- spatial model and directional settings;
- Spencer-taper setting;
- spline knots and orders;
- noise model and spatial weighting;
- `dgsigma` freeze iteration;
- Tukey settings;
- measurement-variance and monotonic settings;
- damping value; and
- convergence and solver diagnostics.

The effective configuration is stored on `dgmodel`. It includes the spatial
model, spline breakpoints and orders, noise and spatial weights,
`dgsigma` freeze iteration, Tukey setting, measurement-variance configuration,
monotonic configuration, Spencer taper, dimensionless damping, the actual
normal-matrix damping scale, convergence tolerance, iteration count, and
convergence status. Solver outcome and error information remain on
`dgsolver_flag`.

Focused tests confirm that the configuration and solver flags survive an
xarray NetCDF save/reload round trip. This concerns fuvpy output content, not
the later icBuilder restart mechanism.

## 5. Run the complete-orbit pilot and select damping — pending (next)

Process 3--5 complete representative orbits with the intended production
configuration. Short extracted sequences are insufficient because temporal
splines and WIC reflattening depend on the complete sequence.

For each orbit:

- inspect `dgmodel`, `dgimg`, `dgweight`, and `dgsigma`;
- inspect iteration counts, convergence, and solver flags;
- check for monotonic fallbacks and artifacts;
- save and reload the output;
- record runtime and approximate memory use; and
- compare a small, defensible range of dimensionless `dampingVal` values.

The outcome is one selected damping value and evidence that the chosen model
behaves acceptably on complete orbits.

## 6. Freeze the fuvpy production configuration — pending

After the pilot, write down the selected settings and damping value, add a
focused regression test for that configuration, and update the public
documentation. At that point fuvpy is ready for downstream icBuilder
integration and archive reprocessing.
