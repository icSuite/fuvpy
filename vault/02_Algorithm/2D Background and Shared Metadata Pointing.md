# 2D Background and Shared Metadata Pointing

Last reviewed: 2026-10-19
Status: background methods and corpus metadata guards tested locally; server rerun pending.
Repository: background_model_asymmetry at46f31e9.

## Public methods and compatibility

The user explicitly requests retaining backgroundmodel_BS with identical1D
defaults, alongside a separate backgroundmodel_BS2D. The new method is exported
at package level and through the xarray fuv accessor. The original1D method's
four fitted fields and numerical diagnostics are bitwise equal to a frozen
synthetic reference generated from committed65603e0. Earlier optional sigma,
prior and initialization controls remain available, disabled by default.

BS2D defaults: cubic x/sine surface,13 x controls, four sine controls,
quadratic120-minute temporal target, numerical damping0.01, Tukey cutoff5,
relative amplitude-change tolerance0.01, variance weighting and quality>=2.
The second axis is exactly solar_view_sine_factor, with its viewing/solar
amplitudes; no raw azimuth, normalization or extra Spencer term is fitted.
Monotonic x increments are constrained independently at every sine/time
control. Tensor coefficients use a sqrt(number_of_sine_controls) basis scale.

## Temporal and independent results

dgsigma_freeze_iteration accepts positive integers below max_iterations,
default2. Pooled residual sigma is estimated after every solve through that
iteration, then the pixel field freezes. At least one subsequent weighted
solve consumes the final estimate before convergence is accepted. Tukey
weights continue updating. The empirical law is sigma²=a*background²+b.

Independent fits run automatically after temporal convergence. Each refits
original counts with fresh all-ones robust weights and the final orbit[a,b]
law evaluated at its own current background. No orbit shape, coefficients,
initial mask or frame-wise noise reestimation is inherited. Sigma is residual
spread, distinct from measurement variance and background-fit uncertainty.

Temporal products retain dgmodel/dgimg/dgweight/dgsigma. Independent products
append _independent. Both coefficient tensors, iteration/convergence flags,
final changes, sigma-law history and penalties are saved. Old background
products are removed from the input copy before fitting, preventing stale1D
decompositions or rerun dimension conflicts. Failure/nonconvergence or all-zero
Tukey support raises; input datasets remain unchanged.

For continuity with the reviewed trials, a coarse legacy base plus tapered
sine reference design sets only the numerical ridge scalar. It is never
fitted/subtracted in BS2D. Default temporal reference targets are120/60minutes;
the2D surface itself uses one120-minute time basis. Per-frame reference uses
a singleton time basis. This reference scale is recomputed from supplied data.

## Shared metadata correction

Both background methods accept a sensor mapping and pointing_correction=True.
Correction requires all three WIC/SI12/SI13 datasets together and occurs once
before per-sensor fitting. DefaultFalse preserves old behavior. Loader option
pointing_metadata=True retains original state, camera angles, exact millisecond
pointing_time and raw native GLAT/GLON/SZA/DZA, and enables viewing factors.
Public correct_pointing can perform this shared step before hemisphere/science
selection; complete-orbit metadata should precede science-frame filtering.

WIC supplies one common correction timeline. Isolated position outliers use
an eleven-frame robust leave-centre-out cubic, threshold max(20km,median+
10MAD-sigma), and local isolation test. Default method=clock_phase also smooths
unwrap(phase)+360*SPIN with a local quadratic; method=position omits smoothing.
Every sensor receives the same orbit solution evaluated at its own precise
time, associated with the nearest WIC timestamp independently of SPIN.
Equidistant matches are left unchanged. Camera alignments stay fixed. WIC
counter transitions other than+1 and time gaps>1.5median cadence retain the
existing segmentation rule. Segments shorter than the eleven-frame window,
gaps and unsupported reference endpoints are not extrapolated. These breaks
affect pointing smoothing only, not the temporal background spline.

Native rays are recovered from saved coordinates and observer metadata. A
full3D Earth-frame rotation is constrained by historical SUN3, stored SZA and
observer vectors. Longitude-only rotation fails baseline parity. This supports
native128x128 SI grids without interpreting their retained256-grid increments.
Changed rays regenerate geographic coordinates, SZA/DZA and solar-view factors.
The packaged original epoch1997 Apex lookup reproduces FUVVIEW3 magnetic
coordinates; MLT retains the original per-frame subsolar magnetic offset.
Baseline reconstruction and magnetic-grid checks reject inconsistent inputs.

Brightness and detector variance stay unchanged. Science coordinates retain
original support intersected with corrected visibility/DZA support. A frame
with no science pixels still receives the shared state repair when raw rays
exist; no invalid pixels are restored. Original state, coordinate hashes,
correction vectors, support and baseline errors are retained. Reapplication
is rejected. This is metadata consistency repair, not independent attitude
validation or an image/limb-driven registration.

## Corpus metadata cases: server orbits197 and365

The user ran icBuilder's read-only diagnostic on the full server corpus.
SI12/197 `s1220002032117.idl` and `s1220002032119.idl` have zero valid stored
and reconstructed DZA pixels. The original strict loader raises for either
frame before returning a sequence, preventing fitting all three sensors.
The user approved excluding these files while retaining the rest of the orbit.
New read_idl option skip_unmapped=True skips only files with neither stored
nor reconstructed valid DZA. It warns with each filename and records a
newline-separated skipped_unmapped_files attribute; calibration arrays remain
aligned. Partial/inconsistent geometry still raises. If every file is
excluded, a clear error reports that no frames remain. DefaultFalse preserves
strict1D loading; icBuilder's2D native writer explicitly enables this option.

WIC/197 repeats22503 at08:25/08:27,22536 at08:28/08:29, and22806 at
21:19/21:20/21:21. WIC and SI13/365 repeat54609 at12:25/12:27. Times increase
and there is no evidence from counter repetition alone of duplicate images.
At21:19:33 on2000-07-21, WIC has22806 while SI12 has22913,48ms apart; the
SI12 frame is one of the unmapped records, so its metadata are not independently
validated. These records nevertheless invalidate an unconditional common
unique-SPIN assumption. The former corrector both validated strict increase
and looked up reference frames by SPIN. The user approved retaining repeated
counter frames and time-based camera association. correct_pointing now
requires finite integer counters without strict increase, finds the nearest
precise WIC time, and evaluates corrections only inside that supported WIC
segment. Original counters and all input frames remain. Existing counter
segmentation and clock-aware phase construction are unchanged.

138 package tests pass, including offset/repeated SI counters, repeated,
reset and jumping WIC counters, equality to separately fitted segments,
unchanged unsupported frames, equidistant timestamp rejection, strict loader
compatibility, individual/all-unmapped cases and retained geometry errors.
Server raw197/365 files are not local; actual corrected products still require
a server rerun. icBuilder's [[IMAGE IDL Pointing Metadata]] owns exact files.

Complete local85 and968 metadata were compared against the committed
pre-fix corrector. Every data variable is bitwise unchanged for all six
sensor/orbit streams (236/226/235 and328/326/326 input frames). The known
85 and968 outliers remain repaired in all cameras. Ignored report:
figures/validation/pointing_counter_time_20261019/complete_orbit_comparison.json.
The native2D writer also completed full85 with all70 saved data variables
per sensor exactly equal to its prior output,122/197/195 science frames and
status1 for all three cameras. icBuilder owns that pipeline evidence.

## Verification and integration boundary

138 package tests pass, including frozen1D defaults, freeze1/2/3 consumption,
known-truth surfaces and abrupt backgrounds with aurora, native SI128 shapes,
joint correction in both methods, input preservation, failure paths and NetCDF.
The built wheel includes both modules and checksum-verified magnetic data.

Complete85 public2D fits converge for122WIC/197SI12/195SI13 frames; SI results
match the saved private trials within1.5e-6counts. WIC's old trial retained a
ridge from earlier paired support; recomputation changes it0.392% and the
background by<1count. With its archived ridge the new temporal solver matches
the saved trial. This is a scalar-reference difference, not a changed surface.

Complete968 joint metadata correction plus public2D fits converge for
312WIC/310SI12/310SI13 frames, with both stages saved/round-tripped. The shared
SPIN25522 position error is repaired in all sensors at their own times.
Complete85 metadata correction also passes all sensors, including WIC's
fully science-masked SPIN3722 while SI retains valid pixels.

Ignored local evidence: figures/validation/2d_shared_pointing_20261019 contains
six native comparison files, real_orbit_regression.json,
joint_pointing_fit_regression.json, implementation_checks.json and the wheel.
The original implementation checkpoint did not change icBuilder entry points.
icBuilder has since integrated this API into an original-style native writer,
retaining its CLI/output/status behavior. The reported server metadata cases
now have local regression coverage; actual197/365 processing remains to
verify on the server. Downstream Product1 field selection remains separate.
