# Production Readiness Audit

Last reviewed: 2026-09-12
Branch: `background_model_asymmetry`

## Conclusion

Do not start the complete IMAGE reprocessing yet. No new numerical defect was
confirmed in the selected directional BS path, and the full local suite has
100 passing tests. The remaining blockers are integration, validation, and
long-run safety.

## Blockers

1. `icBuilder/scripts/pipeline/make_orbit_nc_files.py` still calls the old
   loader and BS configuration. It does not request measurement variance,
   viewing geometry, the directional model, Spencer taper, binned `dgsigma`,
   bin-population weighting, iteration-0 freezing, variance weighting,
   monotonicity, or diagnostics.
2. Production dimensionless damping values have not been selected for either
   BS or SH. Historical numeric values do not retain their old meaning after
   scaling regularization by the initial weighted normal matrix.
3. The selected model has not been evaluated on complete representative
   orbits through the full BS -> SH -> NetCDF round trip. Short example
   sequences do not reproduce the temporal spline or WIC reflattening of a
   complete orbit.
4. The icBuilder runner writes directly to final NetCDF paths, does not verify
   existing outputs for restart, and records the availability array only after
   the complete sensor run. An interruption can therefore leave an invalid
   final file and incomplete run status.

## Repairs before the pilot

- Restore `img_variance` metadata removed from `read_idl`: `excludes`,
  `calculated_from`, and `reflattening_assumption`. The current validation
  script otherwise rejects WIC variance produced by the current loader.
- Enable and inspect `dgfit_iterations`, `dgfit_converged`, and
  `dgsolver_flag` in production outputs.
- Choose a monotonic-failure policy. If ordinary fallback is enabled, correct
  the `monotonic_constraint` metadata when the saved result is no longer
  monotonic.
- Record the actual BS and SH configuration in each output file.
- Add SH convergence diagnostics rather than relying only on a printed
  warning at the iteration limit.

## Known non-blockers

- The unused asymmetric-coordinate helper still overwrites its calculated
  knots with a hard-coded list. It does not affect the selected
  `legacy_x_azimuth` model, but keeps the suite from passing cleanly.
- `monotonic=True` constrains the base component only. The combined base plus
  sine model is not globally monotonic. This is a known model choice, not a
  hidden implementation defect.
- Revised residual norms and replacement of normal equations can wait while
  those diagnostics are unused.

## Verification

- Full suite: 100 passed, 3 failed. The failures are the unused asymmetric
  knot override and two missing measurement-variance metadata fields.
- Directional NetCDF round-trip coverage passes.
- A focused scalar NetCDF check confirmed that the icBuilder compression
  dictionary can save and reload integer, boolean, and floating-point fields.
- Packaged mission-correction tables cover every day from 2000-01-01 through
  2003-07-31 for WIC, SI12, and SI13.

## Next action

Repair the metadata and production runner, then process 3--5 complete
representative orbits through BS, SH, save, and reload. Inspect convergence,
solver flags, scientific fields, runtime, and memory; select BS and SH damping
before launching the complete corpus.
