# FUVVIEW3 calibration tables

These are the original IDL save files used by the recovered FUVVIEW3 source:

- `fuv_mission_correction.sav`
- `wic_flatfield_dbase.idl`
- `si12_flatfield_dbase.idl`
- `si_flatfield_dbase.idl` (SI13)

`fuvpy.src.measurement_uncertainty` reads these files to reconstruct the
mission-correction times flat-field multiplier applied before the IMAGE FUV
files were saved. The implementation follows `fuv_mission_correction.pro` and
`get_flatfield.pro`; the tables themselves have not been altered.
