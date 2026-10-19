# Current State

Last reviewed: 2026-10-19
Branch: `background_model_asymmetry`
Status: original1D preserved; native2D integrated; corpus pointing fixes tested locally

## Latest: separate2D and joint metadata correction implemented

User confirms implementation of a second2D method while preserving existing1D
behavior/defaults. Public backgroundmodel_BS2D automatically converges temporal
x/sine fits then independent frame fits, sharing the configurable final sigma
relation. Temporal dgmodel/dgimg/dgweight/dgsigma remain; independent products
append _independent. Both stages and diagnostics are saved. Sigma-freeze
iteration is any positive integer below max_iterations, default2.

Both methods accept WIC/SI12/SI13 mappings with optional shared metadata
pointing correction, defaultoff. Loader pointing_metadata=True retains state
and exact times. Joint correction repairs one WIC-referenced spacecraft
solution across all cameras and recalculates geographic/magnetic/background
geometry on native grids. No image fitting; counts/variance/alignments retained.

127 tests pass, including bitwise frozen1D defaults and joint1D/2D integration.
Full85 model comparisons and full968 shared-correction+2D fits pass all sensors;
metadata85 also passes. Wheel/source/data packaging checks pass. Details and
local evidence are owned by [[2D Background and Shared Metadata Pointing]].
icBuilder now uses these methods in a separate original-style native writer.
The full server corpus exposes metadata cases absent from local85/86/261/968:
SI12/197 has two entirely unmapped DZA frames; WIC/197 and WIC/SI13/365 repeat
SPIN at different increasing times. WIC and an unmapped SI12 record at almost
the same time also have different counters. At the audit loading/validation
rejected whole sequences, and correct_pointing used SPIN for reference association,
segmentation and phase. The user approved skipping the two unmapped SI12
files, retaining repeated-counter frames and matching cameras by time.
Implemented skip_unmapped=False loader opt-in (enabled by icBuilder's2D
writer), per-file warnings/recorded exclusions, integer-only SPIN validation
and nearest precise WIC time association. Existing WIC counter/time-break
segmentation is retained; short segments keep original pointing and remain
eligible for background fitting. Partial/inconsistent geometry still raises.
138 package tests pass, including counter/segment and loader regressions;
the original1D numerical fixture still passes. Actual197/365 products require
a server rerun because those raw files remain remote. The algorithm note owns
exact evidence and limits.

Complete85/968 metadata results are bitwise identical to the committed
corrector across all three cameras, including both shared position outliers.
The native2D full85 pipeline also exactly preserves all70 saved data fields
per sensor. Reports are linked from the algorithm note;197/365 remain to
verify on the server.

## Earlier opt-in: reestimate sigma after the second fit

`dgsigma_freeze_iteration=2` estimates pooled sigma after fit1 and again after
the Tukey-weighted fit2, then freezes that second sigma field. A following
weighted solve is mandatory even with loose stopping tolerance. Default1 is
unchanged; supplied sigma coefficients cannot be reestimated. Stored[a,b]
describe the last estimate. Tukey weights continue updating normally.
75 focused background tests pass, including default preservation, weighted
reestimation, frozen field, following-solve use and invalid-option checks.

icBuilder's full312-frame orbit calibration exactly reproduces the original
first estimate and verifies final sigma frozen at estimate2. Sigma at modeled
background1000 decreases665→355. Its repaired292 independent-frame comparison
is visually modest, in the same direction but smaller than the earlier manual
quiet-scale trial. Geometry/knots/variance/population weights/ridge remain
fixed; only[a,b] are shared, not an orbit shape or Tukey weights. No production
default selected. icBuilder owns the private frame evidence and proposed
separate SZA/DZA surface with limited azimuth terms; that spatial change is
not implemented.

## Earlier icBuilder direction: initialize from orbit Tukey weights

icBuilder's user rejects the background-shape prior and accepts the known
frame292 metadata position repair. New `initial_robust_weights` optionally
seeds the first single-frame solve with supplied objective weights and sigma.
Every fitted pixel is reweighted afterwards, including initial zeros. The
all-ones damping scale is retained, isolating initialization from persistent
regularization changes. Cannot combine this with an active shape penalty.
67 focused background tests pass, including unchanged all-ones behavior,
recovery of rejected quiet data and mostly-auroral nightside synthetic truth.

The icBuilder fourteen-frame trial with repaired292 geometry finds little
final real-image improvement: reweighting largely removes the initial
benefit. All42 fits converge; initial zeros demonstrably recover. Warm start
alone is not a successful fix for the observed nightside oversubtraction.
Further auroral-rejection strategy is unresolved; no production change.

## Earlier opt-in orbit-background shape constraint (user rejected for icBuilder)

`backgroundmodel_BS` now accepts `background_prior` and nonnegative
`background_prior_strength` for a single frame with supplied empirical sigma.
It adds a fixed population/measurement-weighted prediction-departure penalty
to the robust pixel objective; lambda1 matches the un-rejected data quadratic.
Numerical damping is centered on a feasible projected prior and fresh Tukey
weights initialize against that reference. Monotonic projection honors the
usual fallback flags. Fully rejected data can yield a prior-based prediction.
Zero strength preserves existing outputs exactly. This is empirical shape
regularization, not a calibrated prior variance or covariance calculation.

`store_spatial_coefficients=True` saves per-frame base/directional amplitudes
in `dgmodel_spatial_coefficients` without changing fitted fields. This permits
evaluating the frozen fitted function after a separately estimated geometry
repair. 59 focused background tests pass, including unchanged defaults,
prediction reconstruction, weighted penalty identity, narrow cubic support,
failure handling and moderate broad-aurora attenuation versus strength.

icBuilder owns the fourteen-frame 0968 comparison and separate repaired-position
frame292 trial. All60 independent fits converge; the complete stored baseline
reproduces exactly. Strong regularization restores old background artifacts;
moderate widespread aurora can still be absorbed. No production strength or
uniform improvement/auroral-preservation claim is adopted.

## Opt-in shared residual-scale experiment

For the icBuilder independent-frame background trial, `backgroundmodel_BS`
now accepts `store_sigma_model=True` to save the empirical [a,b] relation
sigma²=a*background²+b in `dgmodel.attrs['sigma_model_coefficients']`.
`sigma_model_coefficients=[a,b]` supplies that relation instead of estimating
frame-bin noise statistics, evaluating it at the current model each IRLS
iteration. Tukey weights initialize afresh. The shared scale remains separate
from measurement variance and background-fit uncertainty.

Default calls retain their existing numerical behavior and output fields.
The focused tests include exact stored/default equivalence, supplied
scale use, fresh weighting, coefficient validation, NetCDF serialization and
known-truth abrupt-background/aurora recovery and monotonic behavior. Eight
icBuilder 0968 preview fits converge; their mixed visual results are owned by
icBuilder's `Background Subtraction - Independent Frame Refinement` note.
No default algorithm, temporal settings or production choice changed.

## Active background architecture

The original background-removal entry point is `backgroundmodel_BS()`; the
separate `backgroundmodel_BS2D()` is described above. The original supports:

- the ordinary one-dimensional BS background; and
- the selected directional model `B0(x) + T(SZA)*S*Bs(x)` through
  `spatial_model="legacy_x_azimuth"`; and
- the experimental coestimated spherical-harmonic correction through
  `spatial_model="legacy_sh"`.

The directional model uses `x=cos(SZA)/cos(DZA)`, the pole-safe solar-view
sine factor `S`, and the optional Spencer taper `T`. The cosine factor is
still calculated and saved for visualization. `monotonic=True` constrains the
base component, not the complete base-plus-directional model.

`backgroundmodel_SH()`, `backgroundmodel_SBS()`, their helpers and public
interfaces, and the rejected compressed-nightside `asymmetric_legacy`
experiment have been removed. New corrected-image workflows use `dgimg`.
Historical NetCDF files with SH fields remain readable as ordinary xarray
datasets.

The new `legacy_sh` option does not restore the old sequential SH solver. It
coestimates a monotonic one-dimensional base and a Spencer-tapered,
non-axisymmetric real SH component in one fit. The historical `n-m` even term
selection is retained while `m=0` is excluded to reduce overlap with the base.
The complete validation-orbit comparison confirms that all twelve WIC/SI fits
are numerically viable with variance weighting after the reviewed problematic
frames are removed. The SH correction often looks competitive and sometimes
better, but substantial base/SH cancellation can still occur and some frames
retain a stronger large-scale gradient. It remains an experimental leading
alternative rather than a selected production model; see
`vault/03_Validation/Coestimated SH Asymmetry Test - 2026-09-14.md`.

## Candidate processing configuration

The current candidate combines:

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

The candidate numerical settings were selected on 2026-09-18: 120-minute base
spacing, 60-minute directional spacing, dimensionless damping `1e-2`,
measurement-variance weighting, and ordinary-solver fallback after a failed
monotonic solve. A controlled orbit-0085 WIC test shows that full-orbit reflattening followed by a
20-frame fit retains the good short-sequence behavior. The poor complete-orbit
result was subsequently isolated to WIC measurement-variance weighting rather
than full-orbit reflattening or the five temporal breakpoints: disabling
variance weighting restores an acceptable long fit. The saved variance remains
useful downstream, but WIC inverse-variance fitting is not production-ready.
A follow-up fit restricted to 00:48--04:51 (122 frames) behaved well with
variance weighting enabled and stayed close to the historical background in
representative early, middle, and late frames. The failure therefore depends
on fitting the heterogeneous complete orbit, not on variance weighting being
intrinsically unusable for every WIC interval.

The final orbit-0085 frames after about 04:49 are saturated garbage images.
They must be rejected before fitting; their large residual contribution can
distort the temporally coupled model even under inverse-variance weighting.
Temporal complexity is now specified by elapsed-time spacing rather than a
fixed knot count. `temporal_knot_spacing` sets the target spacing for `B0` in
minutes; the implementation uses the smallest number of equal intervals no
longer than that target, with both sequence endpoints retained. Directional
models accept `directional_temporal_knot_spacing`; when omitted it inherits
the `B0` spacing. The base and directional terms use separate temporal bases,
so this is a real difference in model freedom rather than metadata only.

The complete-orbit memory blocker is resolved. The solver now exploits the
separable spatial and temporal basis and accumulates `GTG` and `GTd` frame by
frame instead of materializing a whole-orbit pixels-by-coefficients matrix.
The complete 328-frame orbit-0968 WIC asymmetric variance-weighted fit now
finishes in 32 seconds with a measured peak RSS of 6.2 GiB. The previous dense
design matrix alone would have required 13.45 GiB, before several similarly
large working copies.

## Retained revisions

- `dgsigma` is estimated from one statistic per populated geometry bin rather
  than repeating bin statistics per pixel. The obsolete legacy model has been
  removed.
- Every populated geometry bin receives equal total squared influence within
  a frame. The obsolete Gaussian spatial-weighting alternative has been
  removed.
- `dgsigma` is always estimated after the initial fit and then held fixed while
  later model and Tukey updates continue. This is no longer configurable.
- Tukey values enter as objective weights through `sqrt(Tukey)` row
  multipliers.
- BS damping is scaled by the median supported positive diagonal of the
  initial weighted normal matrix.
- `read_idl(..., measurement_variance=True)` attaches detector measurement
  variance before optional WIC reflattening.
- `read_idl()` always assigns `frame_quality`: 0 rejected, 1 detector-valid,
  and 2 science-ready. The initial audited thresholds reject absent mapping
  or relative signal spread below 0.1 and require mapped absolute MLAT to
  reach 68 degrees for science-ready status. Thresholds can be overridden.
- WIC reflattening estimates its shared background from frames with
  `frame_quality >= 1` but applies the correction to every loaded frame.
- `backgroundmodel_BS()` defaults to retaining `frame_quality >= 2` before
  constructing its temporal basis. `minimum_frame_quality=None` disables the
  filter; datasets without the variable retain backward-compatible behavior.
- `variance_weighted=True` applies measurement precision to the BS inversion
  without including it in the empirical `dgsigma` calculation.
- The low-cost monotonic solver has explicit success/fallback flags, and every
  result stores scalar convergence diagnostics.

## Repository cleanup

The obsolete article examples, notebooks, paper figures, generated comparison
figures, Binder environment, empty submodule file, unused `image_orbit.h5`, and
generated caches were removed. Generated figures are now ignored.

The retained research support is:

- `scripts/compare_models.py`;
- `scripts/extract_legacy_background.py`;
- `scripts/validate_background_revision.py`;
- the WIC, SI12, and SI13 IDL validation sequences under `data/`; and
- compact historical BS and SH comparison fields under
  `data/legacy_background/`: `dgmodel`, `dgimg`, `dgweight`, `dgsigma`,
  `shmodel`, `shimg`, and `shweight`.

## Verification

The full suite reports 50 passes. It was consolidated from 88 by combining
repeated error/parameter cases and removing tests of validation-script plotting
plumbing; package-level scientific behavior remains covered. `img_variance` now records what it was
calculated from, which uncertainty sources it excludes, and the deterministic
WIC reflattening assumption. No active source, test, or script refers to SH,
SBS, or `asymmetric_legacy`.

Monotonic solver failures no longer terminate processing. A successful
ordinary fallback is flagged separately from an all-NaN failed product, and
the saved monotonic metadata describes the result actually returned.

Every BS result now records its effective fit configuration and scalar
convergence state on `dgmodel`; solver outcome remains in `dgsolver_flag`.
Focused tests confirm that these values survive NetCDF save/reload.
`backgroundmodel_BS()` always fits `img`; the obsolete alternative `inImg`
path has been removed. It always returns a copy and rejects the obsolete
`inplace`, `diagnostics`, `input_variance`, `monotonic_solver`, `dzalim`, and
`minlat` arguments. Variance-weighted fits always read `img_variance`, and
production monotonic fits always use the low-cost normal-system solver.

The bundled comparison now also plots the coestimated SH model, base, SH
component, corrected image, Tukey weights, and `dgsigma` for all 58 frames
from each sensor.

The complete-orbit SH comparison used the reviewed clean intervals and
measurement-variance weighting. It fitted 122/122/122 WIC/SI12/SI13 frames for
orbit 0085, 190/190/190 for orbit 0086, 272/271/272 for orbit 0261, and
285/283/283 for orbit 0968. All fits succeeded and 2,602 frame panels were
written under `figures/model_comparison_long_sh_variance_filtered/`.

The accumulated normal system agrees with an explicitly constructed dense
design to numerical precision. A complete old-versus-new directional,
variance-weighted monotonic fit agrees to about `2.3e-11` in the modeled and
corrected images. The complete orbit-0261 comparison also finishes for WIC,
SI12, and SI13 and writes all requested variance/no-variance diagnostics.
The same three-way diagnostic analysis was completed for orbit 0968 over the
requested 2001-10-21 14:57 through 2001-10-22 00:38 interval: 284 WIC frames
and 282 frames for each SI sensor.

The frame-quality audit now includes variation within individual detector
rows. `row_relative_spread` stores the coefficient of variation for every
mapped row and `median_row_relative_spread` summarizes each frame. The default
minimum is 0.06; rows require at least eight mapped pixels. Orbit-0968 WIC at
14:00:41 has a median row spread of 0.030,
while the minimum among the other 893 frames accepted by the preceding rules
is 0.083. The new criterion therefore removes the detector-wide row-gradient
frame without discarding another previously accepted validation frame.

Detector-state screening is also part of frame quality. `read_idl()` retains
`hv_phos` and `hv_mcp`; WIC requires both to be at least 1000 V, while SI
requires `hv_mcp >= 1000 V`. Correlation is used only after an observed
detector off-to-on transition. Operational successors remain rejected while
their common-mapped-pixel correlation with the preceding detector-off frame
exceeds 0.90. The common mapped area must cover at least 10% of the detector.
Other quality failures do not start correlation chains.

For orbit 0085, the turn-on rule rejects WIC 00:42--00:46 and makes 00:48:02
the first science-ready frame. With the default loader geometry, correlations
with the voltage-off reference are 0.982, 0.980, and 0.515 at 00:44, 00:46,
and 00:48. Orbit-0968 14:00:41 is rejected independently by its 0.030 median
row spread, while 14:02:43 is science-ready. Turn-on correlations are stored
as `turn_on_reference_correlation`; no support field is stored because the
10% support requirement is fixed metadata on that diagnostic.

## Production boundary

Do not start complete IMAGE processing yet:

1. Complete points 5--6 of
   `vault/03_Validation/FUVPY Production Reprocessing Checklist.md`.
2. Update icBuilder to call the revised loader and BS-only model and to use
   `dgimg`/`dgweight` instead of SH fields.
3. Make the orbit runner interruption-safe with atomic saving and output
   validation.

Detailed decisions, algorithms, and audit evidence remain in the dedicated
vault notes rather than this live status file.

## Complete-orbit comparison data

The external validation corpus is published as Zenodo record 22828032. Its
stable layout separates full-resolution inputs as
`data/idl_files/orbit_####/{wic,si12,si13}_data/` from historical products in
`data/legacy_nc_files/`. The repository ignores the complete `data/` tree.
`scripts/download_data.py` downloads the record's ZIP archive, verifies its
Zenodo MD5 checksum, installs the data in the repository, and removes the ZIP.
`scripts/compare_models.py` checks for both required data subdirectories and
reports the download command when they are absent.

Orbits 0085, 0086, 0261, and 0968 are complete relative to the historical files
for all three sensors.
This corpus remains external because the compressed archive is approximately
5 GB and is unsuitable for Git.

The legacy NetCDF files store dates rather than literal source filenames; their
minute timestamps map directly to the IDL filenames. In orbit 0085 WIC, 88 of
235 saved frames have no pixels in the historical support mask and an entirely
NaN `dgmodel`; they did not contribute rows to the normal equations. The
original code nevertheless passed all 235 northern frames, so their timestamps
set the temporal-knot range and all frames participated in WIC reflattening.
The long script now reproduces that behavior for every available IDL file and
leaves pixel exclusion to `backgroundmodel_BS()`. All 235 historical
orbit-0085 WIC timestamps are now available.
