# Production Readiness Audit

Last reviewed: 2026-09-13
Branch: `background_model_asymmetry`

## Conclusion

Do not start the complete IMAGE reprocessing yet. No new numerical defect was
confirmed in the selected directional BS path, and the consolidated local
suite has 38 passing tests. The remaining blockers are integration,
validation, and long-run safety.

## Blockers

1. `icBuilder/scripts/pipeline/make_orbit_nc_files.py` still calls the old
   loader and BS configuration. It does not request measurement variance,
   viewing geometry, the directional model, Spencer taper, variance weighting,
   or monotonicity.
2. Production dimensionless BS damping has not been selected. Historical
   numeric values do not retain their old meaning after
   scaling regularization by the initial weighted normal matrix.
3. The selected model has not been evaluated on complete representative
   orbits through the full BS -> NetCDF round trip. Short example
   sequences do not reproduce the temporal spline or WIC reflattening of a
   complete orbit.
4. The icBuilder runner writes directly to final NetCDF paths, does not verify
   existing outputs for restart, and records the availability array only after
   the complete sensor run. An interruption can therefore leave an invalid
   final file and incomplete run status.

## Repairs before the pilot

- Inspect the always-saved `dgfit_iterations`, `dgfit_converged`, and
  `dgsolver_flag` fields in production outputs.
- The monotonic-failure policy is implemented. Ordinary fallback is optional;
  all failure paths return explicit flags without terminating processing, and
  `monotonic_constraint` reflects the saved result.
- The actual BS configuration is now recorded on `dgmodel` and verified across
  a NetCDF save/reload round trip.

## Known non-blockers

- `monotonic=True` constrains the base component only. The combined base plus
  sine model is not globally monotonic. This is a known model choice, not a
  hidden implementation defect.
- Revised residual norms and replacement of normal equations can wait while
  those diagnostics are unused.

## Verification

- Current suite: 38 passed after consolidating repeated cases and removing
  validation-script plotting tests.
- Directional NetCDF round-trip coverage passes.
- A focused scalar NetCDF check confirmed that the icBuilder compression
  dictionary can save and reload integer, boolean, and floating-point fields.
- Packaged mission-correction tables cover every day from 2000-01-01 through
  2003-07-31 for WIC, SI12, and SI13.

## Next action

Process 3--5 complete representative orbits through BS, save, and reload.
Inspect convergence, solver flags, scientific fields, runtime, and memory;
select BS damping before updating the downstream production runner.
