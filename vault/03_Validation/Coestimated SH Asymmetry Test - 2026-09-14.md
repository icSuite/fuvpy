# Coestimated SH Asymmetry Test

Last reviewed: 2026-09-15
Branch: `background_model_asymmetry`

## Purpose

Test whether the selected one-dimensional background can be augmented with a
more flexible dawn--dusk correction based on the spherical-harmonic framework
formerly used by the sequential SH correction.

## Implemented model

The optional `spatial_model="legacy_sh"` fits

```text
B = B0(x, t) + T(SZA) H(glat, glon - subsolar_longitude, t)
```

in one inverse problem. `B0` retains the legacy
`x = cos(SZA)/cos(DZA)` spline. `H` is a real, non-axisymmetric spherical-
harmonic expansion. Degree and order default to four. The historical `n-m`
even selection is retained, but `m=0` is excluded because it would overlap
with the coestimated base. The default terms are `(1,1)`, `(2,2)`, `(3,1)`,
`(3,3)`, `(4,2)`, and `(4,4)`, each with cosine and sine components.

The SH columns use the same temporal basis, measurement-variance weights,
bin-population weights, Tukey iterations, and damping as the base. When
requested, monotonicity constrains only `B0`. The combined model and its
components are saved as `dgmodel`, `dgmodel_base`, and `dgmodel_sh`.

## Verification

- A synthetic degree-one test recovers an exact base-plus-SH background and
  confirms `dgmodel = dgmodel_base + dgmodel_sh`.
- SH configuration and missing-coordinate validation are tested.
- The SH output and metadata survive NetCDF save/reload.
- The complete test suite reports 40 passes.
- `scripts/compare_models.py` regenerated all 58 WIC, 58 SI12, and 58 SI13
  bundled example figures with the joint SH model included.

## Preliminary observation

The fit is numerically stable in the inspected WIC, SI12, and SI13 examples,
and the auroral Tukey rejection remains recognizable. The SH component does
capture broad directional structure. However, the plots also show substantial
tradeoff between the SH component and the one-dimensional base even after
excluding `m=0`. The combined corrected image can change little while the two
components individually become large. This is an identifiability warning and
must be evaluated across complete orbits before this model is considered for
production.

## Complete-orbit comparison

`scripts/compare_models_long.py` was extended to compare the coestimated SH
model on the four complete validation orbits with measurement-variance
weighting enabled. WIC is reflattened from the available complete northern
sequence before the reviewed useful interval is selected. Frames with no
BS-supported pixels are then removed entirely from the temporal fit.

The fitted clean-frame counts were:

| Orbit | WIC | SI12 | SI13 |
| --- | ---: | ---: | ---: |
| 0085 | 122 | 122 | 122 |
| 0086 | 190 | 190 | 190 |
| 0261 | 272 | 271 | 272 |
| 0968 | 285 | 283 | 283 |

All twelve orbit/sensor fits completed without monotonic-solver warnings or
failed products. The exhaustive comparison produced 2,602 diagnostic panels
under `figures/model_comparison_long_sh_variance_filtered/`.

Inspection of representative early, middle, and late frames shows that the SH
model remains smooth and temporally usable on complete sequences. In some
frames it represents broad asymmetry more coherently than the retained
sine-spline model. In others the combined corrections are nearly identical,
and the SH base/directional cancellation remains obvious. The earliest retained
orbit-0968 WIC frames also show a stronger residual large-scale gradient under
the SH model, so the result is not uniformly better.

## Status

Experimental leading alternative. It is numerically viable on the complete
validation orbits, but its scientific advantage over `legacy_x_azimuth` is
not uniform and the component cancellation remains unresolved.
