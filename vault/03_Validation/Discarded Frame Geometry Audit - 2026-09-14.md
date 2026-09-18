# Discarded Frame Geometry Audit

Last reviewed: 2026-09-14
Branch: `background_model_asymmetry`

## Question

Do the frames excluded from the orbit-0085, orbit-0086, and orbit-0968
comparison intervals share a viewing-geometry or spacecraft-position
signature? Orbit 0261, for which the complete sequence was retained, provides
a useful comparison.

## Evidence

The audit used WIC fields from the historical orbit NetCDF files and spacecraft
metadata from the corresponding IDL files. The directly relevant metrics are
usable mapped-pixel fraction, maximum visible MLAT, coverage and DZA within
60--80 degrees MLAT, spacecraft geocentric radius and geographic latitude,
mapped-image signal statistics, and the FUVVIEW3 calibration flag.

There are 870 retained and 178 discarded WIC frames across the four validation
orbits. Among the discarded frames:

- 87 have no usable mapped WIC pixels after the loader mask;
- 35 more have mapped pixels but no coverage of 60--80 degrees MLAT;
- 28 more contain that latitude band but its median DZA exceeds 45 degrees;
- 28 pass those coarse geometry checks and require signal- or orbit-specific
  interpretation.

The discarded frames have a median spacecraft radius of 29,741 km, compared
with 47,846 km for retained frames. Their median auroral-band DZA is 44.7
degrees, compared with 31.1 degrees for retained frames. These distributions
overlap, so neither number is a universal threshold.

FUVVIEW3 calibration flag 1 is a strong correlated marker: it occurs in 131 of
178 discarded frames (73.6%) but only 27 of 870 retained frames (3.1%). It is
not a quality flag: FUVVIEW3 defines flag 1 as corrected counts (flat-field and
temperature/sensitivity corrections applied, without conversion to Rayleighs).
The association is therefore empirical, not a sufficient rejection rule or a
causal explanation, and exceptions occur in both directions.

## Orbit-specific findings

- Orbit 0085: the leading WIC sequence is mostly unmapped/zero and flag 1.
  Immediately before the retained interval, the mapped detector signal is
  nearly uniform around 300 rather than containing useful auroral structure.
  Later discarded frames either lie after a large data gap or lose useful
  northern coverage. The exact 04:51 boundary is not a sharp geometric change.
- Orbit 0086: all 23 discarded leading frames occur at geocentric radii of
  9,106--19,379 km. Maximum visible MLAT rises only from about 5 to 68 degrees;
  the first 16 frames do not contain the 60--80 degree band at all. The retained
  sequence begins as auroral coverage and image contrast emerge.
- Orbit 0261: the entire sequence retains polar coverage even though radius
  falls to 21,820 km. This demonstrates that radius alone is not defensible as
  a quality criterion.
- Orbit 0968: the leading discarded frames retain the pole, but the spacecraft
  rises from only 31 to 44 degrees geographic latitude and median auroral DZA
  improves gradually from about 49 to 39 degrees. The 90th-percentile auroral
  DZA remains about 68--64 degrees. There is no physical discontinuity at the
  chosen 14:57 boundary. At 00:41, however, pointing changes abruptly and
  maximum visible MLAT falls from 89.6 to about 33 degrees; all auroral-band
  coverage disappears.

## Conclusion

Most discarded frames lie outside the useful polar-observing portion of the
orbit, but there is no single spacecraft radius, latitude, timestamp position,
or calibration flag that separates all good and bad frames. The mapped
footprint is the correct primary evidence because it already combines position,
attitude, detector field of view, and coordinate mapping.

A future frame-quality rule should separately test:

1. presence and amount of usable mapped data;
2. coverage of the scientific MLAT region;
3. DZA distribution within that region, preferably a high percentile rather
   than only the median;
4. degenerate, saturated, or implausibly uniform detector signal.

The FUVVIEW3 calibration flag should be retained as provenance but not used as
a quality criterion.

These criteria identify whether a frame is useful for auroral science. A
separate decision is needed about whether non-auroral frames should contribute
to the temporally coupled background fit.

## Composite-rule check

The manually selected comparison intervals are not reliable binary truth. A
simple union of objective failure modes was tested without using the calibration
flag:

1. no finite mapped image pixels;
2. maximum mapped MLAT below 68 degrees; or
3. mapped-image coefficient of variation below 0.1, identifying nearly uniform
   frames with little useful structure.

This union identifies 150 of the 178 manually discarded frames. It also rejects
the completely unmapped orbit-0085 frame at 02:17:42, even though that frame lies
inside the manually retained time interval.

The remaining 28 manually discarded frames are the gradually improving leading
geometry of orbit 0968. Rejecting leading frames until the 90th
percentile DZA in the 60--80 degree MLAT band first falls below 64 degrees places
the boundary at 14:57:54 and captures those 28 frames. That exact threshold is
descriptive of the chosen boundary, not independently justified, and should not
be promoted to a production rule without validation on more orbits.

Orbit 0085 at 04:51:09 and orbit 0968 at 00:38:50 were reclassified as retained.
Both have useful mapping, polar coverage, ordinary image contrast, and acceptable
geometry; their earlier exclusion came only from the requested time boundaries.

The defensible interpretation is therefore to use hard per-frame rejection for
missing mapping, insufficient latitude coverage, and degenerate signal, while
treating oblique viewing geometry as a graded quality flag. The old interval
labels should be revised before they are used to validate a production filter.

## Generated diagnostics

- `figures/discarded_frame_geometry_diagnostics.png`
- `figures/discarded_frame_geometry_metrics.csv`
- `figures/frame_quality_audit/orbit_####_wic_frame_quality.pdf`
- `figures/frame_quality_audit/frame_quality_labels.csv`
