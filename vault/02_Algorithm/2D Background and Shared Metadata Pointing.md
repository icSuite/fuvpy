# 2D Background and Shared Metadata Pointing

Last reviewed: 2026-10-19
Status: implemented and verified; icBuilder pipeline integration pending.
Repository: background_model_asymmetry at65603e0 plus uncommitted changes.

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
time and matching spin; camera alignments stay fixed. Missing spins, short
segments, gaps and unsupported reference endpoints are not extrapolated.

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

## Verification and integration boundary

127 package tests pass, including frozen1D defaults, freeze1/2/3 consumption,
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
No icBuilder production entry point or external corpus was overwritten.
Next: integrate the selected outputs into icBuilder after reviewing its field
selection, pointing-before-quality order and dual-result serialization.
