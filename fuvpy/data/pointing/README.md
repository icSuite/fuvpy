# Historical IMAGE/FUV magnetic mapping

`mlatlon.1997a.xdr` is the original lookup table distributed with FUVVIEW3,
copied unchanged from its `support/` directory. FUVVIEW3's `geo_to_apex.pro`
describes this as an Art Richmond Apex-coordinate table valid for epoch1997.
It contains two XDR big-endian float32 arrays, magnetic latitude then magnetic
longitude, each on a361-longitude by181-latitude one-degree geographic grid.

`fuvpy.src.pointing` reproduces that producer's bilinear latitude and circular
longitude interpolation. This preserves the historical coordinate convention;
it is not a modern time-dependent magnetic-field calculation. Magnetic local
time retains the original frame's subsolar magnetic-longitude offset. The
unchanged input coordinates are checked against the table before a repair.

Provenance: recovered local FUVVIEW3 distribution at
`/home/bing/Downloads/fuview3/fuview3/support/mlatlon.1997a.xdr`, 2026-10-19.

SHA256: `54894228da8b60dc21cc5337cf112c4f9e40c26c65d308cf31034562f4e5b549`
