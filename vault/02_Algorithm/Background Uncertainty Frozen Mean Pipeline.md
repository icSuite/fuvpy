# Background Uncertainty Frozen Mean Pipeline

Last reviewed: 2026-09-02
Status: Step 1 complete

## Purpose

This note freezes the deterministic mean-processing configuration for the
background-uncertainty go/no-go experiment. Every later baseline and ensemble
member must use this configuration. Only the fitted fuvpy background may vary.

The freeze is based on source snapshots and file checksums because the fuvpy
uncertainty work and parts of the icBuilder research workspace are not yet
committed.

## Scope of the screening baseline

The executable baseline is the current **bin-first** modular pipeline. This is
appropriate for an initial materiality screen because it can be run without
changing production code. It does not implement the accepted detector-first
publication architecture.

If the background contribution is close to the go/no-go threshold, the
downstream test must be repeated after the detector-first pipeline exists.

## 1. fuvpy mean background processing

Use `scripts/pipeline/make_orbit_nc_files.py` from the frozen icBuilder
snapshot with these settings:

- read the saved FUVVIEW3 IDL files with `dzalim=75` degrees;
- retain northern-hemisphere frames;
- reflatten WIC (`reflat=True`);
- do not reflatten SI12 or SI13 (`reflat=False`);
- fit the B-spline background with spatial knots
  `[-3.5, -0.25, 0, 0.25, 1.5, 3.5]`, `stop=0.01`,
  `n_tKnots=5`, `tukeyVal=5`, and `dampingVal=1e-3`;
- fit the spherical-harmonic residual with `Nsh=4`, `Msh=4`,
  `n_tKnots=5`, `stop=0.01`, `tukeyVal=5`, and
  `dampingVal=1e-4`;
- use WIC `shimg` as the corrected WIC signal; and
- use SI12/SI13 `dgimg` as the corrected SI signals.

The ordinary fuvpy result is calculated with diagnostics and uncertainty
disabled. The experimental uncertainty interface must not change the mean
`dgmodel`, `dgimg`, `shmodel`, `shimg`, robust weights, masks, or residual
spread fields.

## 2. Current Product-1 binning

Use the current Product-1 defaults:

- map every sensor independently to its own regular Cubed-Sphere grid;
- use footprint-overlap binning, not centre assignment;
- use a uniform response over each inferred detector footprint;
- require frame fullness of at least 0.1;
- do not apply the optional line-of-sight brightness correction;
- use WIC `shimg`, SI12 `dgimg`, and SI13 `dgimg`; and
- retain the existing Product-1 count, coverage, within-cell spread, and
  weight definitions unchanged.

The grid geometry is defined in `icbuilder/grids.py`:

- position `(0, 90)`;
- orientation `(0, 1)`;
- nominal extent `50,000 km` by `50,000 km`;
- WIC target resolution `200 km`;
- SI target resolution `400 km`; and
- mapping radius `6501.2 km`, corresponding to 130 km altitude.

The present WIC files have 46 by 46 cells. The source docstring that calls
this the canonical 36-by-36 grid is stale.

## 3. Current SI-to-WIC mapping

For this screening test, use the executable Product-2 mapping in
`icbuilder/precipitationimage.py`:

- direct copy when source and target grids are identical;
- otherwise bilinear interpolation between regular CS-grid cell centres;
- require all four surrounding source values to be finite;
- use nearest-cell inheritance only in the outer half-cell boundary; and
- do not fill internal missing-data gaps.

For standard deviations, propagate the bilinear weights as

\[
\sigma_t = \sqrt{\sum_i a_i^2\sigma_i^2},
\]

assuming independent source-cell errors. This covariance assumption is held
fixed; the go/no-go experiment tests only fitted-background uncertainty.

The area-weighted native-detector SI-to-WIC method in
`icbuilder/detector_coregistration.py` belongs to the accepted detector-first
architecture but is not used by the executable Product-2 pipeline. It is
therefore not silently mixed into this first screening experiment.

## 4. Time support, Kp, proton correction, and precipitation

Keep the following choices fixed:

- match required sensor frames one-to-one within 2 seconds;
- use the latest required-sensor timestamp as the product time;
- obtain definitive three-hour Kp from
  `icbuilder/data/gfz_kp_2000_2003.json` without temporal interpolation;
- use SI12 as the proton-flux source;
- use Hardy et al. (1991) for proton mean energy;
- approximate Hardy corrected geomagnetic coordinates using the Product-2
  Modified Apex grid at 130 km;
- clip proton mean energy to the camera-response range 0.47--46.7 keV;
- apply the existing icPhysics proton response and correction equations; and
- combine method-required sensor weights by multiplication.

Evaluate both precipitation branches without changing their current
equations or masks:

1. `image_ratio`: Frey WIC/SI13 inversion with the existing low-count
   fallbacks (`SI13 < 3`, `WIC < 50`) and E0 clipping;
2. `zhang_paxton`: collapsed Zhang--Paxton E0 with WIC-derived electron
   energy flux.

### Zhang--Paxton execution blocker

The bundled Zhang--Paxton lookup has dimensions `901 x 36 x 36`, whereas the
current WIC Product-1 grid is `46 x 46`. The current Zhang--Paxton branch
therefore cannot run on the frozen Product-1 grid. This pre-existing mismatch
must be resolved before Step 6, either by generating a matching frozen lookup
or by using a verified grid-independent collapsed `(Kp, MLT)` evaluation.

Resolving that mismatch is not part of Step 1 and must not alter the
image-ratio branch or the fitted-background experiment.

## 5. Conductance

Use `icphysics.robinson_conductance` without changing the Robinson equations,
the existing E0--Fe covariance input, or the current zero-flux uncertainty
treatment. Product 3 inherits the precipitation grid, time support, Kp,
proton fields, and combined weight without recalculation.

## Frozen source snapshots

| Repository | Branch | HEAD | State relevant to this test |
|---|---|---|---|
| fuvpy | `background_uncertainty` | `763de63d08676b632a238c67a4e3389c0cb1edad` | Background uncertainty implementation is uncommitted |
| icBuilder | `modular_pipeline` | `aad3f002f5165d642079340e6ab02b618d442ffa` | Large dirty research worktree; current pipeline source files below are frozen by checksum |
| icPhysics | `main` | `c3484fd376caec6616c9476eceab487c7dc718f6` | Clean |
| icReader | `modular_pipeline` | `1fd47caaae8464c9c6d013a329fc409368e56d99` | Clean |
| ZhangPaxton2008 | `main` | `bc0b85be09e123f845b5ab16b0a152125ed0f6f9` | Clean |

### fuvpy checksums

```text
3cf471627fbd5299919a4ad23dc6efeca5da8b25ee4629e0ddd3527e3e7c8d56  fuvpy/src/load.py
9027a0ad401ff3b0bf5c333cd18c1f66b7ff766065e4786c6dac72439ebd7de8  fuvpy/src/background.py
3149ca02fbcf2afa62d1ffd0e4e43ce45182621e58856db9092046ff1bb3c7cb  fuvpy/src/background_fit_outputs.py
34070f513dbf6db3bcbc0253cc2e323764dd11ca9c4b2ee11bbbd2b1e0c7c325  fuvpy/src/background_uncertainty.py
```

### icBuilder checksums

```text
b80fdf8e52c0ee2b5a096fe6b2ab5020a57ff650fbb8f2addc94d11e14fe630b  scripts/pipeline/make_orbit_nc_files.py
3408c52e8e0880ff52c969e87ec013abf3dc4d23c527f079310f3fcaa248e9ac  scripts/pipeline/make_binned_orbit_files.py
d3f2d95e60b4fbd1304be11c4e0eaee1af49b608e21783f301599b3ca41fe2c1  scripts/pipeline/make_precipitation_orbit_files.py
51a54c8cc454870c5dc61f3430f2a151d25ff156e1f9248c4a6feb4f6580bdfa  scripts/pipeline/make_conductance_orbit_files.py
5386ceeae4870fa857d20f303120b6f3fcba45b1e4b1c6dbaf9a8bbf0cf28237  icbuilder/binnedimage.py
18193d34ffea36b6d8515a50d6a6b1eb6014f0c4f70818afb569cb676a85beaf  icbuilder/footprints.py
041fc1e7c0140cb4ef2e8ed13da74220e52982a959bb24e76b3c5b38a9f993e9  icbuilder/precipitationimage.py
68bd1fb492d5e460ee9bb3d819a11e7020926ffe12bd1116b427e81b8e6244d4  icbuilder/grids.py
6c4168ddcd2e510893a4fdaedf70f025db9eadf032129f63355c0fb0c26a56fc  icbuilder/kp.py
29218cf7ab50d7629c92f3fbd2709554dabe593e34f31d82b64bc75b29b2d362  icbuilder/data/gfz_kp_2000_2003.json
36fa6c1bc53ecdb0feaaf83327bd9f4dc64f6864fdacf425d804cc6e31b4a9c8  icbuilder/conductanceimage.py
```

### Shared physics checksums

```text
095527c0ca0330af6a2ee10636545c337f3a08f266add772b88d1f03a3e1eb4d  src/icphysics/image.py
8b27b374640bfb7af3f1912248054fc908a7171896650570bd50b41b1f2be7e0  src/icphysics/hardy.py
fc51c8ce74925d7ce633b0bb0ba94dac75c557e78cd99381078df45a9a01d55a  src/icphysics/hardy_coefficients.py
ea967756e985cb48dd953317b6fce95fb3594fe5c0238dfd307b5260ed23e7e8  src/icphysics/robinson.py
30510480f6b5909e128a7bff9a130be129f5759bd4200f011380ed1c9dff9cae  src/icphysics/zhang_paxton.py
cda87583219160da0c8bed7022de60d0d5e4c1a04dd16dc5d62b52995739b6ad  src/icphysics/zhang_paxton_lookup.py
53c06e2c2ebb5ac8c605185aab7bba65575630be5339bb0d7e1c789262d07e65  src/icphysics/data/zhang_paxton_e0_lookup.nc
```

## Step-1 outcome

The deterministic screening baseline is frozen. No scientific calculation or
pipeline output was generated in this step.

Two limitations must remain visible:

1. the screening baseline is bin-first, while the accepted publication
   architecture is detector-first; and
2. the Zhang--Paxton lookup is incompatible with the current WIC grid and
   cannot be propagated in Step 6 without a separately verified resolution.

