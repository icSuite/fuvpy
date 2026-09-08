# SI12 Measurement Variance Trace

Last reviewed: 2026-09-03
Status: detector-statistics model and opt-in loader integration implemented

## Result

The saved SI12 `image` field can be related directly to the photon-event
counts for the available calibration-flag-1 files. This provides a physical
measurement-variance model for SI12 without estimating variance from the
fuvpy background-fit residuals.

For a detector pixel with raw event count (N_i), define

\[
a_i = g_t f_i,
\]

where (g_t) is the date-dependent mission correction and (f_i) is the
flat-field multiplier. With airglow subtraction disabled and calibration flag
1, FUVVIEW3 saves

\[
y_i = a_i N_i.
\]

The SI detector is a photon counter, so

\[
N_i \sim \operatorname{Poisson}(\lambda_i).
\]

The corresponding variance in the corrected-count units fitted by fuvpy is

\[
\operatorname{Var}(y_i\mid\lambda_i)=a_i^2\lambda_i.
\]

Using the observed count as a plug-in estimate gives

\[
\widehat{\operatorname{Var}}(y_i)=a_i^2N_i=a_i y_i.
\]

This is not simply `variance = image`: the flat-field and mission corrections
change the variance by one remaining factor of (a_i).

## FUVVIEW3 trace

The relevant path through the recovered source is:

1. `getudf_var.pro:423-428` reads the SI event image from
   `return_value.si_pixels`, rotates it, and uses `CONGRID` to place it in the
   internal 256-by-256 structure.
2. `getudf_var.pro:500-510` records calibration flag 1, obtains the SI12
   flat-field, and reads the date-dependent mission correction.
3. `getudf_var.pro:533` calculates
   `cal_img = (raw_image * temp_correction) * ff_corr` for the current files,
   because their stored airglow scale is zero and their Rayleigh conversion is
   one.
4. `downshift.pro:63` returns the SI image to 128 by 128 before the IDL save.

The mission correction is defined in `fuv_mission_correction.pro:23-30`. It
normalizes SI12 sensitivity to mission day 171 of 2000 using stellar
measurements. The flat-field selection and multiplication convention are in
`get_flatfield.pro:108-247`. Before and after the 2000 boom loss are treated
as separate periods. After the boom loss the row profile is also rotated and
shifted, so the effective multiplier can be two-dimensional.

## Direct check using the saved files

The 40 SI12 example frames from 2000 day 137 all have:

- `CALIBRATION_FLAG = 1`, meaning corrected counts;
- `AIRGLOW_SCALE = [0, 0, 0]`; and
- a five-second SI viewing exposure according to the instrument paper.

For that date:

- mission multiplier: 0.990879;
- usable-pixel combined multiplier (a_i): 0.990879--1.500695;
- median combined multiplier: 1.178129.

Dividing the saved values by (a_i) recovered integer event counts in all
532,623 positive usable samples. The largest distance from an integer was
(4.7\times10^{-6}) count. This also demonstrates that the intermediate
128-to-256-to-128 operations do not average independent event counts in the
final usable pixels for this pre-boom sequence.

The reconstructed total count rate over usable pixels was approximately
22,100--25,700 counts per second, below the published SI system limit of
100,000 counts per second. This does not exclude local count-rate effects in
the brightest structures, but it argues against global detector saturation in
this example.

## Literature constraints

Mende et al. describe SI12 and SI13 as crossed-delay-line single-photon
counters. Each detected photoelectron is assigned to the TDI memory image.
They give:

- 128-by-128 memory pixels;
- a five-second viewing exposure per 120-second spacecraft rotation;
- expected SI12 sensitivity of 0.018 count per Rayleigh per pixel per
  exposure;
- maximum system rate of 100,000 counts per second;
- detector-background rate below 5 counts per second over the photosensitive
  area; and
- response uniformity of approximately plus or minus 20 percent.

The background rate corresponds to fewer than 25 detector events across a
whole five-second image, or roughly 0.0015 event per 128-by-128 pixel if
uniform. It is already part of the measured Poisson count. The response
uniformity is handled by the flat-field and is a calibration systematic; it
must not be added as independent pixel measurement noise.

Frey et al. report an in-flight SI12 sensitivity of 23.5 corrected counts per
kR for a five-second exposure. This is consistent with the interpretation of
the saved calibration-flag-1 values as sensitivity-normalized counts.

## Remaining implementation choices

The physical model is sufficiently specified for positive, unsaturated SI12
pixels. Before connecting it to the background covariance calculation, three
small choices remain:

1. Port the FUVVIEW3 SI12 flat-field selection, post-boom rotation, and mission
   correction into a readable Python helper, or preserve the effective
   multiplier when IDL files are converted.
2. Choose a nonzero estimate of the Poisson mean for observed zero-count
   pixels. A small count floor or a Poisson posterior mean avoids assigning
   exact zero variance. Only about one percent of usable pixels were zero in
   the checked sequence.
3. Mask pixels where detector saturation or invalid calibration factors make
   the linear Poisson transformation unreliable.

The resulting covariance must be labelled detector-count or measurement
statistics conditional on the supplied calibration. It excludes uncertainty
in the flat-field, mission-gain correction, absolute calibration, and the
background model itself.

## Evidence paths

- FUVVIEW3 source: `/home/bing/Downloads/fuview3/fuview3/`
- SI12 instrument paper:
  `/home/bing/Dropbox/work/code/repos/icBuilder/literature/Mende_et_al_2000_IMAGE_SI.pdf`
- quantitative calibration paper:
  `/home/bing/Dropbox/work/code/repos/icBuilder/literature/Frey_et_al_2003_IMAGE_FUV_quantitative_interpretation.pdf`
- checked SI12 files:
  `/home/bing/Dropbox/work/code/repos/icBuilder/example_data/s12_data/`
