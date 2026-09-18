# 2D Model Experiments - Summary

Last reviewed: 2026-09-18
Source branch: `2D_model` at `1556ff1`
Status: durable findings transferred; branch is an experimental archive

## Question

The legacy coordinate

```text
x = cos(SZA)/cos(DZA)
```

cannot distinguish viewing directions with the same SZA and DZA. The branch
tested whether retaining more viewing geometry could represent dawn-dusk
background asymmetry without fitting aurora.

## Models tested

The experiments covered:

- free two-dimensional SZA-DZA spline surfaces;
- signed DZA and eastward viewing-angle coordinates;
- pole-safe solar-view cosine and sine factors;
- separate base, cosine, and sine spline surfaces;
- alternative illumination coordinates, monotonic constraints, staged robust
  fitting, and several terminator tapers; and
- reduced directional models using SZA, the legacy coordinate, or a single
  scalar asymmetry coefficient.

## Durable findings

- Reconstructed DZA matched the stored FUVVIEW3 DZA to about `1.5e-5` degrees
  on the tested WIC, SI12, and SI13 frames. The viewing-geometry reconstruction
  itself is therefore sound.
- Unsigned DZA cannot distinguish dawn from dusk. Directly signing DZA creates
  a seam, while an eastward-angle coordinate becomes ill-defined at the pole.
- Pole-safe solar-view dot and oriented-cross-product factors preserve the
  relative solar/viewing direction without an angular wrap or polar
  singularity. Both remain useful diagnostic fields.
- Free two-dimensional and multi-component surfaces were too flexible. Their
  components traded against one another, absorbed auroral structure, and
  produced ringing or unstable robust weights.
- Constraining only a base component did not make the combined directional
  background monotonic. The cosine component was particularly prone to
  unwanted structure and brought no consistent visual benefit.
- Alternative illumination coordinates, nightside compression, and more
  elaborate sine-coordinate surfaces did not consistently improve corrected
  images enough to replace the legacy base coordinate.
- Smoothly suppressing directional freedom toward the nightside was useful.
  The normalized Spencer spherical-Chapman profile was the most defensible
  taper tested, although it remains optional.

## Retained formulation

The strongest candidate was the low-freedom sine-only model

```text
B = B0(x) + Tspencer(SZA) S Bs(x),
x = cos(SZA)/cos(DZA),
directional_x_knots = [-3.5, 0, 1.5, 3.5].
```

It keeps the established one-dimensional base, adds one signed directional
profile, and avoids a free SZA-DZA surface. This formulation, the optional
Spencer taper, and both solar-view diagnostic factors were reimplemented
cleanly on `background_model_asymmetry`; the larger model zoo was not carried
forward.

## Evidence boundary

The branch included geometry checks and focused synthetic tests, but model
selection was driven mainly by visual inspection of selected frames rather
than representative complete-orbit validation. The experiments identify a
promising constrained formulation, not proof of physical correctness or final
knot tuning. Reopen the discarded 2D variants only if new validation evidence
shows a specific limitation of the retained model.
