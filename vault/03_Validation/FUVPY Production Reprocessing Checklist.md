# FUVPY Production Reprocessing Checklist

Last reviewed: 2026-09-18
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

All results store `dgfit_iterations` and `dgfit_converged`; all-NaN failure
products record `dgfit_converged=False`. `dgsolver_flag` stores the solver
outcome and retains the solver error message in its metadata.

## 4. Make each BS result self-describing — complete

Save the configuration needed to understand and reproduce each background fit,
including at least:

- spatial model and directional settings;
- Spencer-taper setting;
- spline knots and orders;
- fixed binned noise model and bin-population spatial weighting;
- fixed initial-fit `dgsigma` estimation;
- Tukey settings;
- measurement-variance and monotonic settings;
- damping value; and
- convergence and solver diagnostics.

The effective configuration is stored on `dgmodel`. It includes the spatial
model, spline breakpoints and orders, the fixed noise and spatial weights,
initial-fit `dgsigma` estimation, Tukey setting, measurement-variance configuration,
monotonic configuration, Spencer taper, dimensionless damping, the actual
normal-matrix damping scale, convergence tolerance, iteration count, and
convergence status. Solver outcome and error information remain on
`dgsolver_flag`.

Focused tests confirm that the configuration and solver flags survive an
xarray NetCDF save/reload round trip. This concerns fuvpy output content, not
the later icBuilder restart mechanism.

The production interface has also been narrowed: BS always returns a copy,
variance weighting always uses `img_variance`, and monotonic fitting always
uses the low-cost normal-system solver. Experimental or obsolete selectors are
no longer accepted.

Frame screening is now integrated without using FUVVIEW3 calibration state.
`read_idl()` stores one ordinal `frame_quality`: 0 rejected, 1 detector-valid,
and 2 science-ready. WIC reflattening estimates its scalar background from
quality-1-or-higher frames, while BS defaults to quality 2. The current
defaults reproduce the audited hard criteria (mapping present, relative
spread at least 0.1, and mapped absolute MLAT reaching 68 degrees). Their
archive-wide behavior remains part of the pilot rather than a frozen final
selection.

## 5. Complete the orbit pilot and select production settings — in progress

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

Preliminary orbit-0085 WIC diagnostic (2026-09-13): loading and reflattening
all 235 historical frames, then fitting only the same 20-frame interval used
by `scripts/compare_models.py`, recovers the good short-sequence result. Against
the bundled comparison reflattening, the corrected-image median difference is
below one count and the Tukey-weight RMS difference is 0.0035. The poor
complete-orbit revised fit is therefore associated with the temporal fit over
the full duration, not with full-orbit WIC reflattening. Five temporal
breakpoints span about 9.5 minutes per interval in the short test but about
three hours per interval over the complete orbit. Temporal resolution must be
tested before damping is selected.

A subsequent full-orbit ablation identified WIC measurement-variance weighting
as the immediate cause of the unacceptable revised fit: setting
`variance_weighted=False` restores an acceptable result with the same temporal
knots. The reference WIC variance is proportional to the observed count, so
inverse-variance fitting favors darker pixels and frames. In orbit 0085 the
median supported-pixel count is about 4754 in the first half and 6216 in the
second; after bin-population weighting, the first half receives about 1.55
times more inverse-variance influence per supported pixel. Visible off-Earth
stars are excluded by the geometric support mask and are probably a marker of
the darker viewing geometry rather than the cause. Do not enable WIC variance
weighting in production until observed-count weighting has been justified or
replaced. The variance field can still be saved for downstream propagation.

A second ablation retained variance weighting but restricted the orbit-0085
fit to 00:48--04:51. The 122 available frames span 00:48:02--04:49:10. Both
the revised base and asymmetric fits completed, and the revised background
remained close to the historical background in representative early, middle,
and late frames. This narrows the problem to the interaction between inverse
variance weighting and the heterogeneous complete-orbit fit. It does not
support rejecting variance weighting solely from a short or moderately long
interval test.

Inspection identified saturated garbage frames after about 04:49 as the
material difference between the successful restricted fit and the failing
complete-orbit fit. Add explicit saturation quality control before drawing a
production conclusion about variance weighting. Elapsed-time temporal spacing
is now implemented: `temporal_knot_spacing` controls `B0`, while the optional
`directional_temporal_knot_spacing` controls `BS` and otherwise inherits the
base spacing. The code retains both endpoints and divides the duration into
equal intervals no longer than the requested target. Production spacing still
requires a complete-orbit comparison after frame screening.

Frame quality now includes a minimum median within-row relative spread of
0.06. This rejects orbit-0968 WIC 14:00:41, whose apparent whole-image
variation is almost entirely a top-to-bottom row gradient. Correlation is used
only after an observed detector off-to-on transition: operational successors
remain rejected while their correlation with the preceding detector-off frame
exceeds 0.90. This keeps the orbit-0085 00:42--00:46 turn-on rejection without
changing another accepted orbit-0085 frame. Orbit-0968 14:00 is rejected by
the independent row metric and 14:02 is science-ready without a correlation
check. Both row-spread diagnostics and `turn_on_reference_correlation` are
saved for archive-wide review. Rows require at least eight mapped pixels before
their spread enters the frame median, and a turn-on correlation requires at
least 10% common mapped detector coverage.

The historical orbit-0085 BS field is effectively monotonic in the legacy
spatial coordinate to within 1% of its model range in 125 of 147 fitted frames.
The remaining 22 frames, concentrated near the end of the sequence, have
large non-monotonic excursions. The revised coefficient-level constraint is
therefore inactive-looking over most individual frames but still removes real
freedom used by the historical full-orbit solution. A full-orbit ablation of
monotonicity must precede adding temporal knots.

Complete-orbit memory validation (2026-09-14): the BS solver now accumulates
the normal system frame by frame using the separable spatial and temporal
bases. The 328-frame orbit-0968 WIC asymmetric variance-weighted fit completes
in 32 seconds with 6.2 GiB peak RSS. Dense-versus-accumulated algebra and full
model outputs agree to numerical precision. Complete orbit-0261 comparisons
also finish for WIC, SI12, and SI13. This removes the memory blocker but does
not complete point 5: saturation screening, temporal-spacing selection,
damping selection, and scientific inspection remain.

Orbit-0968 interval comparison (2026-09-14): the full three-way diagnostic
set was generated from 2001-10-21 14:57 through 2001-10-22 00:38 for WIC,
SI12, and SI13. All requested variance, no-variance, and direct-difference
figures completed without solver or memory failure. This is validation input
for point 5, not yet a scientific selection of the production settings.

### Remaining production decisions

The selected spatial model remains `legacy_x_azimuth`; `legacy_sh` is an
experimental comparison and is not a production blocker. The numerical
production candidate was selected on 2026-09-18:

```python
backgroundmodel_BS(
    images,
    spatial_model="legacy_x_azimuth",
    directional_x_knots=[-3.5, 0, 1.5, 3.5],
    temporal_knot_spacing=120,
    directional_temporal_knot_spacing=60,
    spencer_taper=True,
    variance_weighted=True,
    monotonic=True,
    monotonic_fallback=True,
    dampingVal=1e-2
)
```

An archive-wide preflight of the existing MLAT and relative-spread thresholds
was considered but is not required before the run. Before archive-wide
processing, complete the following:

1. Run the exact selected configuration once on the four validation orbits,
   inspect representative and failure-flagged frames, and verify NetCDF
   save/reload before freezing the settings.

## 6. Freeze the fuvpy production configuration — pending

After the pilot, write down the selected settings and damping value, add a
focused regression test for that configuration, and update the public
documentation. At that point fuvpy is ready for downstream icBuilder
integration and archive reprocessing.

The existing directional NetCDF round-trip test now uses the selected
120/60-minute spacing, `dampingVal=1e-2`, variance weighting, monotonic fitting,
and fallback configuration. The README has been corrected to remove obsolete
`input_variance`, `monotonic_solver`, and optional-diagnostics interfaces and
now documents the selected candidate and current failure flags. Final freeze
therefore depends only on the exact four-orbit run and review in point 5.
