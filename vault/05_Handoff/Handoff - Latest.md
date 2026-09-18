# Handoff - Latest

Last updated: 2026-09-18
Branch: `background_model_asymmetry`
Checkpoint: frame-quality filtering implemented; sine-spline experiment removed

## Obsolete branch reconciliation

`background_model_revision` is a Git ancestor of this branch. `experimental`
is not an ancestor, but it has no vault files absent from the current branch;
its durable uncertainty work was selectively carried into the revised
implementation and documentation. The unique `2D_model` experiments are
condensed in `vault/03_Validation/2D Model Experiments - Summary.md`. Their
retained geometry, Spencer taper, and sine-only legacy-x model were
reimplemented here; the remaining model zoo is archival and should not be
merged wholesale.

## What is implemented

The maintained background function is `backgroundmodel_BS()` with either the
ordinary one-dimensional model or the selected directional model:

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

The model is `B0(x) + T(SZA)*S*Bs(x)`. It saves the base and directional
components, both solar-view factors, and the Spencer taper. Monotonicity
applies only to the base component.

Temporal freedom is now set in minutes. `temporal_knot_spacing` controls the
base component and defaults to 120 minutes.
`directional_temporal_knot_spacing` optionally controls the directional
component; if omitted, it inherits the base spacing. Break points are evenly
distributed over the actual sequence duration so the resulting interval is
never longer than the requested target. The old `n_tKnots` argument has been
removed. Both requested and actual spacings and both sets of break points are
stored on `dgmodel`.

`read_idl()` now assigns one ordinal `frame_quality` variable to every frame:
0 rejected, 1 detector-valid, and 2 science-ready. The default criteria use
mapped-pixel presence, relative signal spread 0.1, and mapped absolute MLAT
68 degrees; callers may override these thresholds. WIC reflattening uses only
quality-1-or-higher frames to estimate its scalar background, then applies the
correction to all frames. `backgroundmodel_BS()` defaults to
`minimum_frame_quality=2` and filters before temporal knots are constructed.
FUVVIEW3 `calibration_flag` remains separate provenance.

An experimental alternative is now available as
`spatial_model="legacy_sh"`. It jointly fits the same one-dimensional base and
a Spencer-tapered, non-axisymmetric real spherical-harmonic component in
geographic latitude and solar-relative longitude. It is not the removed
sequential SH correction and does not change the selected production model.

SH and SBS implementations, accessors, exports, helpers, and examples are
gone. The rejected compressed-nightside `asymmetric_legacy` path is also gone.
The corrected image produced by active processing is `dgimg`.

## Repository layout after cleanup

The article-specific `examples/` content and tracked generated figures were
removed because they describe an obsolete algorithm. The retained scripts are:

- `scripts/compare_models.py` — compare the historical BS fields, revised
  one-dimensional model, and selected directional model;
- `scripts/extract_legacy_background.py` — rebuild the compact historical
  comparison files from complete old orbit files; and
- `scripts/validate_background_revision.py` — compare cumulative BS revisions
  on complete native-grid orbits.

The IDL validation sequences and compact historical BS/SH comparison fields
remain under `data/`. Each compact file contains `dgmodel`, `dgimg`,
`dgweight`, `dgsigma`, `shmodel`, `shimg`, and `shweight`. Generated figures
are ignored. The obsolete Binder environment,
empty `.gitmodules`, unused 89 MB `image_orbit.h5`, and local build/test caches
were removed.

## Verification

`pytest -q` reports 50 passes. The suite was reduced from 88 by consolidating
repeated cases and removing tests of validation-script plotting plumbing while
retaining package-level scientific and serialization coverage. `img_variance` now records `calculated_from`,
`excludes`, and the WIC `reflattening_assumption`. The metadata fix does not
change the variance values or background fit. No active source, test, or
script refers to SH, SBS, or `asymmetric_legacy`.

`monotonic_fallback` now controls whether an ordinary fit is attempted after
a monotonic failure. Failed fits return flagged all-NaN products rather than
terminating processing. Flags distinguish successful fallback, no fallback,
and failed fallback, and the saved monotonic metadata reflects the returned
fit.

`dgmodel` now stores the effective BS configuration, including spline
breakpoints and orders, weighting and variance settings, monotonic settings,
damping value and scale, iteration count, and convergence state. The metadata
and solver flags survive NetCDF save/reload.

`backgroundmodel_BS()` now always fits `img`. The obsolete `inImg` option and
its redundant metadata were removed.

`dgsigma` is always estimated after iteration 0 and then held fixed while the
Tukey weights continue to update. The obsolete `dgsigma_freeze_iteration`
option was removed.

The binned `dgsigma` model and bin-population spatial weights are now the only
implementations. The obsolete `noise_model` and `spatial_weighting` options,
legacy pixel-expanded noise model, and Gaussian weighting path were removed.

The remaining unused public controls were also removed. BS always returns a
copy, always saves convergence diagnostics, reads `img_variance` when
`variance_weighted=True`, and uses the low-cost normal-system monotonic solver.
The obsolete `inplace`, `diagnostics`, `input_variance`, `monotonic_solver`,
`dzalim`, and `minlat` arguments are rejected. The full monotonic solver remains
an internal numerical reference.

The coestimated SH experiment has synthetic recovery, validation, component
identity, and NetCDF round-trip coverage. `scripts/compare_models.py`
successfully regenerated 58 WIC, 58 SI12, and 58 SI13 example panels. The
inspected fits are stable, but the base and SH components can trade substantial
amplitude while leaving the combined correction similar. Treat this as an
identifiability warning and do not promote the model without complete-orbit
evidence. Details are in
`vault/03_Validation/Coestimated SH Asymmetry Test - 2026-09-14.md`.

The complete-orbit test is now done with variance weighting and the reviewed
problematic frames removed. All twelve orbit/sensor fits succeeded. The clean
frame counts were 122/122/122 for orbit 0085, 190/190/190 for orbit 0086,
272/271/272 for orbit 0261, and 285/283/283 for orbit 0968
(WIC/SI12/SI13). The run generated 2,602 panels under
`figures/model_comparison_long_sh_variance_filtered/`. The SH correction is
often competitive and sometimes visually better, but the improvement is not
uniform and base/SH cancellation remains visible. It should be compared more
systematically before replacing the selected `legacy_x_azimuth` model.

## Downstream incompatibility

icBuilder still calls `backgroundmodel_SH()` and expects `shimg` and
`shweight`. It cannot use this fuvpy branch until its orbit runner and detector
reader are changed to the BS-only fields.

## Next action

Use `vault/03_Validation/FUVPY Production Reprocessing Checklist.md` as the
numbered source of truth. Points 1--4 are complete. Point 5, the complete-orbit
pilot and damping selection, is next.

A controlled orbit-0085 WIC test loaded and reflattened all 235 historical
frames, then selected the same 20-frame interval used by the bundled short
comparison before fitting. This recovered the good short-sequence result;
Tukey weights differed by only 0.0035 RMS and the corrected-image median by
less than one count. Full-orbit reflattening is therefore not the cause of the
poor long fit. Point 5 must first resolve temporal resolution: five breakpoints
span about 9.5 minutes per interval in the short test and about three hours per
interval over the complete orbit.

The next ablation resolved the immediate failure more specifically: the full
orbit looks acceptable when `variance_weighted=False`. WIC variance is
proportional to observed counts, so inverse-variance fitting favors darker
pixels and changes frame influence across the orbit. For orbit 0085, the first
half has about 1.55 times the effective inverse-variance influence per
supported pixel of the second half. Off-Earth stars are outside BS fit support
and are not the direct cause. Preserve `img_variance` for propagation, but do
not make WIC variance weighting a production default without a revised and
validated statistical treatment.

A follow-up orbit-0085 fit retained variance weighting but restricted the
sequence to 00:48--04:51 (122 frames, actual timestamps
00:48:02--04:49:10). Both revised and asymmetric fits completed. Inspection
of representative early, middle, and late frames shows that the revised base
background stays close to the historical saved background. The variance issue
is therefore conditional on the heterogeneous complete-orbit fit rather than
an inevitable failure on every WIC sequence. Diagnostic figures are in
`figures/model_comparison_orbit_0085_0048_0451_variance/wic/`.

The likely cause is now visible: the final orbit-0085 frames after roughly
04:49 are saturated garbage images. Production fitting needs a saturation
mask/frame-rejection rule before the variance-weighting decision is revisited.
The temporal spline now uses elapsed-time spacing. Production validation still
has to compare defensible base and directional spacings after unusable frames
are rejected.

The complete-orbit memory blocker is resolved. The solver no longer builds the
dense whole-orbit design. It evaluates supported spatial rows by frame and
uses spatial/temporal separability to accumulate the small normal system
(`GTG`, `GTd`) directly. Algebraic comparison with the former dense design is
equal to numerical precision, and a complete old-versus-new model comparison
differs by at most about `2.3e-11` in the image fields.

The formerly failing 328-frame orbit-0968 WIC asymmetric variance-weighted fit
now completes in 32 seconds with 6.2 GiB peak RSS. The complete orbit-0261
variance/no-variance analysis also completed for WIC, SI12, and SI13, writing
272, 271, and 272 plots respectively to each of its three comparison groups.
The generated figures are under `figures/model_comparison_orbit_0261_*` and
remain ignored by Git.

Orbit 0968 was subsequently analyzed over the requested interval
2001-10-21 14:57 through 2001-10-22 00:38. Available frames span 14:57:54 to
00:36:48: 284 WIC frames and 282 frames for each SI sensor. The standard
variance, standard no-variance, and direct asymmetric variance-difference
groups are complete under `figures/model_comparison_orbit_0968_1457_0038_*`.
WIC was reflattened using the complete loaded orbit before the interval was
selected for fitting.

A cross-orbit audit of discarded comparison frames is recorded in
`vault/03_Validation/Discarded Frame Geometry Audit - 2026-09-14.md`. Most
discarded frames have no usable WIC mapping, lack auroral-latitude coverage, or
view that region obliquely. No universal spacecraft-radius cutoff exists;
direct mapped-footprint metrics are more defensible. FUVVIEW3 flag 1 is strongly
associated with discarded frames but has exceptions in both directions. Flag 1
means corrected-count calibration, not bad quality, so this association must not
be used as a causal explanation or automatic rejection rule.

A composite check without calibration flags rejects frames with no mapped data,
maximum mapped MLAT below 68 degrees, or mapped-image coefficient of variation
below 0.1. That audit led to a second spatial-variation metric: the coefficient
of variation is calculated within every mapped detector row and its frame
median must be at least 0.06. Orbit-0968 WIC 14:00:41 has a median of 0.030;
the minimum among the other 893 frames accepted by the preceding rules is
0.083. Per-row and median diagnostics are retained in the loaded dataset.

Detector-state validation is required for quality 1 and 2. `read_idl()` stores
`hv_phos` and `hv_mcp`; WIC requires both voltages to be at least 1000 V and SI
requires `hv_mcp >= 1000 V`. Correlation screening is limited to observed
detector off-to-on transitions. Operational successors remain rejected while
their correlation with the preceding voltage-off frame over common mapped
pixels exceeds 0.90; at least 10% of the detector must be common. The single
saved diagnostic is `turn_on_reference_correlation`. This keeps the established
orbit-0085 00:42--00:46 rejection and makes 00:48:02 the first science-ready
frame. Unrelated quality failures do not propagate to later frames. Thus,
orbit-0968 14:00 is rejected independently by its row-spread metric and 14:02
is science-ready without requiring a recovery comparison.

The full-orbit inputs and historical NetCDF products are published as Zenodo
record 22828032. `scripts/download_data.py` retrieves and verifies the roughly
5 GB ZIP, installs it as the ignored repository-root `data/` directory, and
deletes the downloaded archive. `scripts/compare_models.py` now checks for
`data/idl_files` and `data/legacy_nc_files` before processing and directs the
user to the downloader when either is absent.

Orbits 0085, 0086, 0261, and 0968 are complete for all sensors relative to
their historical files. The
apparent extra SI files in orbit 0085 are southern frames and are not among the
legacy northern inversion timestamps.

The WIC orbit-0085 NetCDF contains another 88 dates with all-NaN `dgmodel`.
Those pixels were excluded from the historical normal equations by the SZA
support mask, but the original pipeline still passed all 235 northern frames:
their dates determine the temporal-knot range, and they participate in WIC
reflattening. The long comparison now does the same and lets
`backgroundmodel_BS()` exclude unsupported pixels internally. All 235
orbit-0085 WIC timestamps are available. The historical and revised fits both
retain five temporal knots.

The production checklist was reconciled on 2026-09-18. The selected numerical
candidate uses 120/60-minute base/directional spacing, `dampingVal=1e-2`,
variance weighting, monotonic fitting, and ordinary-solver fallback. Turn-on
frames that still match the detector-off reference use the saved mapped-pixel
0.90 correlation diagnostic; other quality failures remain independent.
The remaining validation step is to run the exact selected configuration across
the four validation orbits once more.

After points 5--6 are complete:

1. Update icBuilder to use the revised loader and BS-only background call.
2. Add atomic saving and restart validation to the orbit runner.

The production-readiness audit is in
`vault/03_Validation/Production Readiness Audit - 2026-09-12.md`.
