# Current State

Last reviewed: 2026-09-08
Repository snapshot: `background_uncertainty` at `763de63`
Worktree state: uncertainty implementation, tests, documentation, and project
memory are uncommitted

## Confirmed implementation

- The BS fit now has an experimental one-dimensional
  `geometry_model="effective_illumination"`. Its coordinate is the logarithm
  of a smooth solar-illumination term multiplied by a spherical van Rhijn
  viewing factor raised to a configurable exponent. Defaults are 130 km
  emission height, viewing exponent 1, and a 5 degree twilight width. Fixed
  knots emphasize the terminator. The legacy coordinate remains the default.
  Effective illumination defaults to neutral spatial weighting because the
  inverse-Gaussian multiplier is specific to the legacy coordinate.
  The fitted coordinate is saved as `dayglow_coordinate` on the detector
  grid, including the coordinate model and parameter values as attributes.
  The coordinate and all default knots are translated so SZA 90 degrees and
  DZA 0 degrees define zero. The unit-effective-illumination knot is therefore
  positive; this translation changes labels but not the fitted model.
- A second experimental option,
  `geometry_model="bounded_effective_illumination"`, retains the same solar
  and viewing geometry but applies a shifted `log1p` transform. It keeps the
  nadir terminator at zero while compressing the deep nightside toward
  `-log(2)`, so the coordinate changes more gradually around zero. The
  original raw-log effective-illumination option is unchanged. Separate
  comparison figures are stored under `figures/IRLS_iterations_bounded/`.
- `geometry_model="linear_effective_illumination"` retains the smooth solar
  term and finite spherical viewing factor without a logarithmic transform.
  The nadir terminator and subsolar point are normalized to zero and one.
  Default nadir knots at SZA 75, 60, 45, and 30 degrees distribute spline
  freedom across the dayside; oblique dayside knots cover DZA 60, 75, and 90
  degrees. All 129 example-frame comparisons completed under
  `figures/IRLS_iterations_linear/`. The earlier coordinate implementations
  and their figures remain available.
- The `dgsigma` calculation now treats `effective_illumination` and
  `linear_effective_illumination` like the legacy geometry: one bin contains
  the full negative-coordinate nightside, followed by 14 equal-width bins
  from zero to the positive coordinate endpoint. This changes only the
  empirical residual-spread estimate; spline knots and population weighting
  are unchanged. Three focused bin/noise tests passed.
- `geometry_model="cos_sza"` is now the parameter-free geometry comparison.
  Its coordinate is exactly `cos(SZA)` and is independent of DZA. Its default
  spline knots span -1 to 1 with one broad nightside interval, and its
  `dgsigma` calculation uses one nightside bin plus 14 equal dayside bins.
  Four focused geometry/bin tests passed, and 129 comparison figures were
  generated under `figures/IRLS_iterations_cos_sza/`.
- `geometry_model="powered_legacy"` generalizes the legacy coordinate to
  `cos(SZA)/cos(DZA)**viewing_exponent`. Exponent zero reproduces `cos_sza`
  and exponent one reproduces the legacy coordinate. Its default knots retain
  the legacy internal values while covering the powered DZA range, and
  its `dgsigma` calculation uses one nightside plus 14 dayside bins. Five
  focused tests passed. The p=0.5 comparison produced 129 figures under
  `figures/IRLS_iterations_powered_legacy_p0p5/`.
- An unverified `geometry_model="powered_asymmetric_legacy"` applies the same
  viewing exponent before compressing the negative coordinate with the
  existing `nightside_width` transformation. It has been implemented and
  documented, but deliberately not executed or tested yet.
- Two unverified spherical legacy options use `cos(SZA)` multiplied by the
  finite spherical-shell viewing factor. `spherical_legacy` leaves both sides
  uncompressed; `spherical_asymmetric_legacy` compresses only the negative
  coordinate with `nightside_width`. Both reproduce the legacy dayside at
  DZA zero and neither uses `viewing_exponent` or `twilight_width`. Code,
  documentation, and unexecuted tests are present; no comparison figures have
  been generated.
- `geometry_model="asymmetric_legacy"` leaves the positive legacy coordinate
  unchanged and maps its negative side through
  `nightside_width*tanh(x/nightside_width)`. The default width is 0.2, so the
  deep nightside saturates near -0.2 without introducing a discontinuity at
  the terminator. The original legacy model remains unchanged. All 129
  example-frame comparisons completed under
  `figures/IRLS_iterations_asymmetric_legacy/`.
- `backgroundmodel_BS(..., monotonic=True)` now provides an experimental
  spatially monotonic mean fit. It represents each spatial control-point
  sequence as one free baseline followed by non-negative increments, applied
  independently to every temporal basis component. The temporal spline is
  unchanged. The ordinary unconstrained fit remains the default.
  `monotonic_solver="normal"` reduces the pixels to a small normal system
  before bounded fitting and is the default monotonic solver;
  `monotonic_solver="full"` retains the more numerically stable full-design
  BVLS calculation.
- Ordinary BS and SH calls retain their legacy numerical fields and avoid
  uncertainty work. Full-array fingerprints protect fitted models, corrected
  images, weights, spread fields, and masks.
- Opt-in diagnostics retain coefficients, basis metadata, convergence, rank,
  conditioning, and the exact input variance used.
- Fixed-weight BS and conditional SH coefficient covariance propagate to
  fitted-background and corrected-image standard deviations. The corrected
  image retains its covariance with the input used to fit the background.
- Coefficient-space samplers preserve spatial and temporal fit correlation
  after an xarray/NetCDF round trip.
- Repeated-fit routines rerun the complete BS fit or the dependent WIC BS-to-SH
  chain, allowing robust weights and `dgsigma` to change.
- Compact temporal fields store covariance between detector-grid mean fitted
  backgrounds in each frame. They do not claim to be dense pixel covariance.
- The interface is package-level but explicitly experimental. It excludes
  calibration, flat-field, atmospheric-response, regularization-bias, and
  background-model discrepancy uncertainty.

## Complete-orbit evidence

Validation used detector-grid orbit 214 as a bright case and orbit 1926 as a
moderate case for WIC and SI13.

- Repeated/fixed fitted-background sigma ratios were generally about
  1.0--1.15; corrected-image ratios remained near one.
- Median fixed one-sigma coverage of repeated corrected images was about 0.70.
- In the count-scale orbit-214 test, analytic/repeated adjacent-frame
  correlations were 0.979/0.978 for SI13 and 0.980/0.984 for WIC.
- The repeated robust background can be shifted relative to the unperturbed
  fit even when its spread is well approximated. Count-scale median shifts
  were -0.095 SI13 counts and +0.730 WIC counts; covariance alone does not
  encode this nonlinear estimator mean shift.
- `dgsigma**2` is not a suitable detector variance by default. WIC `dgsigma`
  was roughly 3--5.3 times `sqrt(img)` at the median and produced much larger
  nonlinear shifts.
- In the count-scale full WIC chain, median BS/SH component sigmas were
  0.238/4.283 counts and the combined sigma was 4.301 counts. Pixelwise
  BS--SH correlation changed sign spatially (median 0.025), while the combined
  frame-mean background had adjacent-frame correlation 0.995.

The fixed-weight covariance is therefore a useful compact conditional spread
and short-lag correlation estimate under a stated variance model. The full
dependent repeated chain remains the reference for WIC. Long-lag correlation
from only 10--20 refits is too noisy to interpret quantitatively.

## Verification

- The effective-illumination coordinate distinguishes DZA at SZA=90 degrees
  and passed a direct fitted-background monotonicity test. A representative
  20-frame variance-weighted WIC fit converged in three IRLS iterations in
  1.498 s; the minimum ordered background step was zero. The focused
  BS/reference suite passed 32 tests, with only the existing NetCDF binary
  compatibility warning.
- The monotonic option converged on the synthetic WIC sequence, its fitted
  background was non-decreasing with `cos(SZA)/cos(DZA)` in every frame, and
  the protected unconstrained reference output remained unchanged. On one
  20-frame WIC sequence, the normal/full solvers took 1.630/7.958 s and their
  maximum fitted-background difference was (1.52\times10^{-8}). The focused
  BS/reference suite passed 31 tests.
- Conditional covariance is deliberately rejected for `monotonic=True`.
  Active inequality constraints make the existing unconstrained linear
  covariance calculation invalid.
- 18 focused tests passed before the final temporal-review corrections. The
  review then added stricter fixed-mask, covariance-amplitude, positive-
  semidefinite, sample-reconstruction, and NetCDF assertions; these final
  edits were not rerun at the user's request.
- Complete-orbit validation figures are in `/tmp/fuvpy_validation/`; they are
  exploratory products and are not tracked.

## Current boundary and next decision

The software and representative complete-orbit validation through the
experimental interface are complete. No icBuilder product has changed.

The production handoff decision is now to save the background-corrected image
and its estimated measurement variance, while treating the selected
background correction as deterministic. Explicit `dgmodel_fit_sigma` or
`shmodel_fit_sigma` will not be propagated into icBuilder. Background-model
discrepancy remains an unquantified systematic that can be assessed through
method sensitivity tests rather than folded into the measurement variance.
The repeated-fit routines remain validation tools inside fuvpy, not part of
the production data path.

The earliest data in scope are the saved FUVVIEW3 IDL images. The project will
not seek Level-1 UDF files or modify FUVVIEW3 to export detector statistics.
The required input is therefore an **effective post-FUVVIEW observation
variance in corrected-count space**, not a claim of recovered raw detector
variance.

The SI12 detector-count path has now been traced through the recovered
FUVVIEW3 code and calibration tables. For calibration-flag-1 files with no
FUVVIEW3 airglow subtraction, the saved corrected count is (y_i=a_iN_i),
where (N_i) is the photon-event count and (a_i) is the mission-gain times
flat-field multiplier. Its physical measurement variance is therefore
(a_i^2\lambda_i), with plug-in estimate (a_i y_i). Integer raw counts were
recovered to floating-point accuracy from 532,623 usable samples in the
complete 2000-day-137 example sequence. See [[SI12 Measurement Variance
Trace]].

This replaces the residual-fit route for SI12. The SI and WIC variance
equations are implemented in `fuvpy/src/measurement_uncertainty.py` and are
now available through the opt-in `read_idl(...,
measurement_variance=True)` path. The loader reconstructs the FUVVIEW3
mission-gain and flat-field multiplier, preserves it as
`calibration_multiplier`, and adds `img_variance` before optional WIC
reflattening. The original calibration tables are packaged unchanged under
`fuvpy/data/calibration/`. Ordinary loading remains unchanged. For WIC, use
the working interpretation that `IMAGE` is
the flag-1-equivalent corrected-AD field and flag 3 labels the additional
background-subtracted `CIMAGE` product. The literature then supplies a
candidate reference noise scale: 23 equivalent output quanta per 100 R per
2-by-2-pixel resolution cell, together with 612.6 AD/kR per detector pixel.
This is enough for a first-order compound-Poisson model, but not an exact
pixel covariance: distributing the resolution-cell noise among pixels and
neglecting additive CCD noise remain explicit approximations. SI13 remains
unresolved. The complete derivation and uncertainty boundary are recorded in
[[IMAGE FUV Measurement Uncertainty - SI and WIC]].

The new loader path passed syntax checks and one-frame SI12 and WIC smoke
checks; the complete test suite was not run. The SI12 result matched
`calibration_multiplier * img` exactly, and the WIC reference result matched
`10.65 * img` exactly when reflattening was disabled. With reflattening, the
variance remains tied to the pre-reflattened detector image under the stated
frozen-background approximation.

The BS fit now has a separate opt-in `variance_weighted=True` path. It requires
`input_variance`, multiplies each fit row by (1/\sqrt{v_i}), and calculates
Tukey weights from either ((y_i-G_i m)/\sqrt{v_i}) or the legacy empirical
`dgsigma` scale. `tukey_weights=False` disables robust reweighting without
disabling spatial or measurement weighting. The legacy fit remains the default.
Non-finite or non-positive variances at fitted pixels are rejected;
the current measurement model should therefore be loaded with a positive
Poisson offset (for example `poisson_offset=0.5`) before weighted fitting.
Conditional covariance uses the same combined row weights as the fit.

The loader now applies that offset to WIC as well as SI. For reference WIC,
`poisson_offset=0.5` adds the equivalent-count floor
\(0.5(10.65)^2\simeq56.7\) AD squared. The default offset remains zero, so
ordinary and earlier opt-in results are unchanged.

`backgroundmodel_BS()` now also accepts `sza_max`. Values above this limit are
clipped only when constructing the geometry coordinate
`cos(SZA)/cos(DZA)`, so the choice affects the spline design and `dgsigma`
without modifying the stored SZA field. The default remains no clipping.

### Legacy `dgsigma` audit

The published method says to regress the binned residual RMSE against binned
background intensity in the model
\(\sigma^2=s_d^2 I_{BS}^2+\sigma_d^2\). The released `_noisemodel()` instead
copies each bin RMSE to every pixel and regresses those repeated values against
the individual pixel model values. This implicitly weights bins by their pixel
population and introduces within-bin pseudo-replication. All three available
WIC example sequences collapsed to an effectively zero slope, producing
constant `dgsigma` values of about 486, 429, and 1849 counts despite strongly
varying background models.

A new opt-in `noise_model="binned"` implementation now uses one weighted
`model**2` mean and one residual variance per bin. `noise_model="legacy"`
remains the default for reproducibility. In the two 20-frame examples,
`dgsigma` spanned 423--743 and 356--681 counts instead of being constant. The
resulting background models differed from legacy by median absolute values of
about 10 and 7 counts, with maxima of about 286 and 220 counts. This supports
new implementation retains the legacy calculation as an explicit comparison
option. The paper explicitly prescribes one nightside bin,
14 dayside bins, nonnegative variance coefficients, and squared row weights;
those elements are not audit findings to change.

A separate discrepancy remains in the spatial coverage weights. The paper and
its figure script use inverse observed population in narrow geometry bins for
each frame. The package currently uses
`exp(0.5 * (fraction / 1.5)**2)`, introduced in a commit whose source comment
called it a normal-distribution test. It is the inverse of an assumed normal
density, not a normal weight: it upweights large absolute geometry coordinates,
and its square enters the normal equations.

The Gaussian multiplier and scale-freezing choice are now exposed without
changing their defaults. The clearer `spatial_weighting` input now accepts
`"gaussian"`, `"bin_population"`, or `"none"`. The bin-population option
calculates populations separately for each frame and gives every populated
coordinate bin equal total least-squares weight. Because the solver squares
row multipliers, each pixel receives a multiplier proportional to the inverse
square root of its bin population. A common normalization preserves the mean
squared multiplier and prevents the fixed damping from overwhelming the fit.
The earlier literal inverse-count implementation drove a representative
effective-coordinate WIC background to a median of 0.84 counts; the corrected
weighting gave 1055 counts, close to 1068 counts without spatial weighting.
The old `gaussian_weights` Boolean remains available for compatibility but
cannot be combined with the new input. `dgsigma_freeze_iteration=0` estimates
the scale after the initial fit and freezes it; the default value 1 reproduces
fuvpy by recalculating it once after the second fit. These settings are
recorded in fit diagnostics.

The current `IRLS_iterations` figures were generated with the misspelled
`gaussian_weight=False`, so the package silently retained inverse-Gaussian
weighting. They still show a repeatable later-iteration Tukey-weight structure
near `cos(SZA)/cos(DZA)=0`, sometimes visible in the corrected image, but must
be regenerated with `spatial_weighting="none"` before separating this behavior
from the spatial-weight choice.

Code inspection narrows the likely mechanism. Under the applied DZA limit,
fraction zero is the physical solar terminator. The scale gradient is expected
because `dgsigma = sqrt(a * dgmodel**2 + b)`, while the sharp feature can arise
when coherent terminator residuals cross the redescending Tukey cutoff and are
fed back through later scale recalculation. Zero is also the asymmetric
boundary between one broad negative-fraction noise bin and the narrow positive
bins, so binning may affect the scale coefficients, but the returned scale law
is continuous and cannot alone create a sharp spatial line. The cubic spline
has a simple knot at zero and remains continuous there.

Candidate empirical scales can now be inspected without changing that frozen
operational fit. With diagnostics enabled,
`dgsigma_diagnostic_iterations=5` stores the first five candidate
`dgsigma` fields and the empirical Tukey weights each candidate would imply.
The iteration history is explicitly diagnostic: candidates after
`dgsigma_freeze_iteration` are not fed back into the inversion. Setting the
new control to zero disables the history. This separates the question of when
the fit should freeze its scale from the question of how the scale would have
evolved. The focused background-uncertainty suite passes (26 tests).

Five focused tests passed for legacy opt-in behavior, both Tukey-scale choices,
the no-Tukey path, weighted covariance, and rejection of zero fit variance.
The comparison script now
compares the legacy and variance-weighted fits, but its full figure set has not
been regenerated. Because measurement precision changes the scale of the data
term, the old numerical damping value is not directly comparable and must be
assessed before scientific use.

A previous frame-by-frame comparison was generated from all 43 example
IDL files for each of WIC, SI12, and SI13. Files were split into three
contiguous orbit segments before fitting. Each ordinary BS fit was compared
with `uncertainty=True` using `input_variance="img_variance"`; the fitted BS
model, corrected image, and robust weights were numerically identical in the
inspected results, as required. WIC used `reflat=False`. The 129 figures are
under `figures/measurement_variance_bs_comparison/` and can be regenerated by
`examples/scripts/compare_measurement_variance_bs.py`. Those existing figures
show the earlier post-fit-only calculation; rerunning the updated script will
replace that comparison with the variance-weighted fit.

This run also exposed and fixed a post-boom calibration edge detail: IDL
`ROT` extrapolates from the nearest input pixel when `MISSING` is omitted. The
Python flat-field rotation now uses nearest-edge extrapolation instead of
inserting zeros.

Before any icBuilder integration, validate the reconstructed multiplier over
complete SI12 and SI13 sequences and test the WIC equivalent-quanta model.
Then propagate a small background ensemble
through representative processing and compare its contribution with existing
count/binning, coregistration, E0, and conductance uncertainty. Keep the new
component only if it materially affects corrected WIC/SI13, image-ratio E0, or
Hall/Pedersen conductance. Background-model discrepancy remains excluded by
the stated assumption that the selected model is correct.

## Go/no-go test progress

Steps 1 and 2 are complete. The executable bin-first mean pipeline, all
relevant parameters, repository snapshots, and source checksums are frozen in
[[Background Uncertainty Frozen Mean Pipeline]]. The accepted detector-first
architecture is not yet runnable and is not being silently substituted into
this screening test. The current 46-by-46 WIC grid is also incompatible with
the bundled 36-by-36 Zhang--Paxton lookup; that mismatch must be resolved
before downstream Step 6.

Step 2 used ten well-separated common orbits and the final fitted fields from
complete-orbit products. Candidate effective variances are
\(19.0424B\) for WIC, \(0.3317+1.3407B\) for SI12, and \(2.5569B\) for SI13,
all in counts squared. These laws are not yet validated. WIC has strong
detector residual correlation (median 0.75 at lag one and 0.46 at lag five),
while SI correlations are weaker but nonzero. WIC also shows large
orbit-to-orbit changes in residual center and spread. See
[[Effective Observation Variance - Step 2 Results]].

Step 3 is complete and the simple background-only variance model failed the
statistical-validity gate. On held-out orbits 480, 800, 1200, and 1450, pooled
standardized spread was 0.611 for WIC, 1.010 for SI12, and 0.645 for SI13.
SI12 passed all screening checks. WIC and SI13 failed both transferability
between orbits and remaining background/SZA trend checks; SI13 also failed the
DZA trend check. See [[Effective Observation Variance - Step 3 Results]].

Under the agreed plan, do not proceed to downstream ensemble generation or
icBuilder integration with these laws. The next decision is whether to close
the empirical effective-variance route as a no-go for WIC and SI13 or develop
instrument-based models. The immediate bounded next step is to implement and
validate the traced SI12 Poisson model before attempting either remaining
sensor.
