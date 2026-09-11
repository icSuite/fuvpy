# 2D Model Experiments - Summary

Last reviewed: 2026-09-11
Experimental branch: `2D_model`
Experimental checkpoint: `1556ff1`

## Why the experiments were run

The legacy coordinate `cos(SZA)/cos(DZA)` collapses illumination and viewing
geometry to one number. It cannot distinguish geometrically different pixels
with the same collapsed coordinate and therefore has limited ability to fit
dawn-dusk background asymmetry. The experimental branch tested whether a
direction-aware background model could improve this without fitting aurora.

## Broad findings

### Viewing geometry

- SZA and unsigned DZA alone cannot distinguish opposite horizontal viewing
  directions.
- Assigning a sign directly to DZA introduced discontinuities when the
  eastward line-of-sight component changed sign.
- An eastward projected angle removed that seam but depended on the local
  east direction and became unsuitable near the geographic pole.
- Pole-safe solar-view dot and oriented-cross-product factors retained the
  relative solar/viewing direction without choosing geographic east at the
  pole.
- Reconstructed DZA agreed with stored FUVVIEW3 DZA to about `1.5e-5` degrees
  on bundled WIC, SI12, and SI13 frames. This validates the geometry
  reconstruction, not the scientific background model.

### Background surfaces

- A free two-dimensional SZA-DZA surface could alter broad asymmetry, but it
  introduced many coefficients and did not distinguish dawn from dusk.
- Adding independent cosine and sine SZA-DZA surfaces supplied that direction
  but was too flexible. Components could trade against one another, absorb
  auroral structure, and produce ringing.
- Constraining only the base surface to be monotonic did not make the combined
  directional model monotonic.
- The cosine component was repeatedly associated with unwanted structure,
  particularly toward the nightside. Removing it was more promising.
- Nonlinear two-dimensional models using the sine factor as a spline
  coordinate produced poor corrected images or unstable weights in the
  inspected examples.
- A single scalar asymmetry was tested as a strongly constrained alternative,
  but it was not selected over the x-dependent directional profile.

### Coordinates and tapering

Several alternative illumination coordinates and knot layouts were tried.
None consistently improved the corrected images enough to replace the legacy
coordinate. The useful change was instead to retain the legacy base and add a
limited signed directional correction.

Hard or raised-cosine suppression near the terminator could create edge
effects. The normalized Spencer spherical-Chapman profile provided a smooth,
physically motivated taper. It was useful for suppressing directional freedom
on the nightside, but its inclusion should remain optional so its impact can
be tested directly.

## Leading candidate

Visual comparison favored

```text
B = B0(x) + Tspencer(SZA) S Bs(x)
```

with

```text
x = cos(SZA)/cos(DZA)
directional_x_knots = [-3.5, 0, 1.5, 3.5]
```

This was exposed experimentally as
`spatial_model="legacy_x_azimuth"`, `include_cosine=False`, and
`spencer_taper=True`. It retains the successful one-dimensional base, adds
only one signed directional spline, and avoids a free SZA-DZA surface.

The cosine factor remains scientifically useful as a visualization and
geometry diagnostic even though it is not part of the retained fit.

## Evidence level

The branch accumulated focused geometry, synthetic-recovery, solver, taper,
and compatibility tests. At the final checkpoint, 55 focused tests covering
the wider model collection passed. Exact synthetic recovery showed that the
retained formulation and taper were assembled as intended.

Most scientific judgments came from visual inspection of bundled or selected
frames. Representative complete-orbit validation was not run. Consequently:

- the directional model is a leading candidate, not a validated default;
- the proposed knot locations are provisional;
- improved appearance does not establish that aurora is excluded from the
  fitted directional component; and
- the test count demonstrates implementation consistency, not physical
  validity.

## What carries forward

The clean branch should carry forward only:

- the validated viewing-geometry reconstruction;
- both solar-view geometry factors;
- the optional Spencer taper;
- the sine-only legacy-x directional model;
- separate base and directional outputs; and
- focused tests for those pieces.

All other model variants remain available in `2D_model` as experimental
history and should not be merged into the revised background model.

