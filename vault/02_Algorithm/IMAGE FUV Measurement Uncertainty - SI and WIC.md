---
title: "IMAGE FUV Measurement Uncertainty for SI and WIC"
subtitle: "Detector-to-saved-image equations and conditional background-fit propagation"
date: "2026-09-03"
geometry:
  - a4paper
  - margin=24mm
fontsize: 10pt
mainfont: Inter
sansfont: Inter
monofont: DejaVu Sans Mono
mathfont: Latin Modern Math
colorlinks: true
linkcolor: blue
urlcolor: blue
header-includes:
  - |
    \usepackage{booktabs}
    \usepackage{array}
    \usepackage{microtype}
    \usepackage{enumitem}
    \setlist{nosep}
    \setlength{\parskip}{0.45em}
    \setlength{\parindent}{0pt}
---

# Purpose and scope

This note defines candidate measurement-variance models for the IMAGE FUV
SI12, SI13, and WIC images supplied to the `fuvpy` background fits. It answers
four questions:

1. What is the detector-level stochastic model?
2. How does FUVVIEW3 transform the detector output into the saved image?
3. Which variance equation should be evaluated in the saved-image units?
4. Which uncertainty terms are retained or deliberately omitted?

The target quantity is the conditional measurement covariance of the saved
detector-grid image. It is not the total uncertainty of an auroral brightness
or conductance product. Calibration uncertainty, an incorrect background
model, atmospheric emission modelling, proton correction, coregistration,
binning, electron-energy retrieval, and conductance conversion are separate
uncertainty components.

The current working interpretation of WIC is:

- `IMAGEINFO.IMAGE` contains the flag-1-equivalent corrected-AD image;
- flag 3 labels the additional background-subtracted `CIMAGE` product; and
- `CIMAGE` is not the image fitted by the present fuvpy pipeline.

This interpretation is supported by the contents of the saved files and by
fuvpy's conversion of `IMAGE` using 612.6 AD units/kR, but the producer-side
flag-3 code is absent from the recovered FUVVIEW3 source. It therefore remains
an explicit assumption.

# Common calibration equation

Let $x_i$ be the uncalibrated detector output in pixel $i$. For calibration
flag 1, the relevant FUVVIEW3 operation is

\[
 y_i = a_i x_i,
 \qquad
 a_i = g_s(t)f_{s,i}(t),
\]

where:

- $y_i$ is the saved corrected-count or corrected-AD value;
- $g_s(t)$ is the mission sensitivity correction for sensor $s$ and date
  $t$; and
- $f_{s,i}(t)$ is the selected flat-field multiplier.

The more general code path is

\[
y_i =
\frac{[(x_i-b_i)g_s(t)]f_{s,i}(t)-c}{r_s},
\]

where $b_i$ is a FUVVIEW3 airglow estimate, $c$ is a scalar background
level, and $r_s$ converts corrected counts to Rayleigh. The files checked for
the present uncertainty work have zero FUVVIEW3 airglow scale, while flag 1
sets $c=0$ and $r_s=1$. The simpler equation $y_i=a_i x_i$ is therefore
the applicable one. Files that do not satisfy these conditions require their
own propagation and must not silently use the equations below.

The correction inputs come from:

| Quantity | Source |
|---|---|
| Sensor and saved image $y_i$ | `IMAGEINFO.INST_ID` and `IMAGEINFO.IMAGE` in each IDL file; exposed as `id` and `img` by `fuvpy.read_idl()` |
| Time $t$ | `IMAGEINFO.TIME`; exposed as the `date` coordinate |
| Calibration mode | `IMAGEINFO.CALIBRATION_FLAG` in the IDL file |
| FUVVIEW3 airglow settings | `IMAGEINFO.AIRGLOW_SCALE` in the IDL file |
| Mission correction $g_s(t)$ | `support/fuv_mission_correction.sav`, selected by `fuv_mission_correction.pro` and normalized to mission day 171 of 2000 |
| Flat field $f_{s,i}(t)$ | `wic_flatfield_dbase.idl`, `si12_flatfield_dbase.idl`, or `si_flatfield_dbase.idl`, selected by `get_flatfield.pro` |

FUVVIEW3 chooses the nearest available flat field while keeping pre- and
post-boom-loss periods separate. The post-boom field is rotated and shifted.
The multiplier is therefore allowed to depend on both date and detector
location.

The saved IDL files contain the calibration flag and airglow scale. The
opt-in `fuvpy.read_idl(..., measurement_variance=True)` path now preserves
these fields, reconstructs $a_i$ from packaged copies of the original
FUVVIEW3 tables, and adds it as `calibration_multiplier`. No UDF data are
required.

# SI12 and SI13

## Exact conditional variance

The SI detectors are crossed-delay-line single-photon counters. Let $N_i$
be the number of events assigned to output pixel $i$ during one five-second
exposure. Under the usual photon-counting model,

\[
N_i\mid\lambda_i \sim \operatorname{Poisson}(\lambda_i).
\]

For a flag-1 image with no FUVVIEW3 background subtraction,

\[
y_i=a_iN_i.
\]

Conditional on the expected event count $\lambda_i$,

\[
\operatorname{E}(y_i\mid\lambda_i)=a_i\lambda_i,
\]

and

\[
\boxed{
\operatorname{Var}(y_i\mid\lambda_i)=a_i^2\lambda_i.
}
\]

The multiplier is squared because a deterministic multiplication by $a_i$
multiplies a variance by $a_i^2$.

## Plug-in estimate from a saved image

The unknown mean count may be estimated by

\[
\widehat{\lambda}_i=N_i=\frac{y_i}{a_i}.
\]

Substitution gives

\[
\boxed{
\widehat{\operatorname{Var}}(y_i)=a_i y_i.
}
\]

Only one factor of $a_i$ remains after replacing $\lambda_i$ by $y_i/a_i$.
Consequently, `variance = image` is incorrect whenever the flat-field or
mission correction differs from one.

For SI12, this relation was checked directly on 40 files from 2000 day 137.
Dividing 532,623 positive saved pixels by $a_i$ recovered integer event
counts with a maximum numerical discrepancy of $4.7\times10^{-6}$ count.
The combined multiplier ranged from 0.990879 to 1.500695, with median
1.178129. The same FUVVIEW3 detector and calibration structure applies to
SI13, but the integer-recovery check has not yet been repeated for SI13.

## Zero counts and masks

The direct plug-in estimate gives zero variance when $y_i=0$. The true
Poisson mean need not be zero merely because zero events were observed. Two
transparent options are:

\[
\widehat{\lambda}_i=\max(N_i,N_{\min}),
\]

or, using the Jeffreys-Poisson posterior mean,

\[
\widehat{\lambda}_i=N_i+\tfrac12.
\]

The latter gives

\[
\widehat{\operatorname{Var}}(y_i)
=a_i y_i+\tfrac12a_i^2.
\]

This choice has not been made. Invalid calibration factors, known bad pixels,
and pixels affected by saturation or count-rate nonlinearity should be masked
rather than assigned a finite variance.

## SI terms retained and omitted

Retained:

- photon-event counting noise;
- the deterministic scaling of that noise by mission and flat-field
  correction; and
- detector background events insofar as they are already counted in $N_i$.

Omitted:

- uncertainty in $g_s(t)$ and $f_{s,i}(t)$;
- absolute radiometric-calibration uncertainty;
- dead-time, saturation, and other count-rate nonlinearities outside masked
  regions;
- any FUVVIEW3 background-estimation uncertainty when the airglow scale is
  nonzero; and
- background-model discrepancy in fuvpy.

The published SI detector background is below five events per second over the
whole photosensitive area. These events are already part of the observed
Poisson count. They should not be added again as an independent variance term.
The quoted detector-response nonuniformity is a calibration systematic handled
by the flat field, not independent measurement noise.

# WIC

## Detector model

WIC is not a direct integer counter. A photon produces a photoelectron, the
microchannel plate amplifies it, the phosphor and fibre taper produce a light
spot, and the CCD converts the spot into AD units. Approximately 300 CCD video
frames are accumulated using time-delay integration during the ten-second
exposure.

Let $N_i\sim\operatorname{Poisson}(\lambda_i)$ be the number of detected
events contributing to one output element. Let $A_{ij}$ be the uncalibrated
AD amplitude contributed by event $j$, and let $\epsilon_i$ contain
additive CCD, readout, dark, and quantization noise. Then

\[
x_i=\sum_{j=1}^{N_i}A_{ij}+\epsilon_i.
\]

Write

\[
\mu_A=\operatorname{E}(A),
\qquad
m_{2,A}=\operatorname{E}(A^2),
\qquad
\kappa=\frac{m_{2,A}}{\mu_A}.
\]

Assuming independent event amplitudes and zero-mean additive noise,

\[
\operatorname{E}(x_i)=\lambda_i\mu_A,
\]

\[
\operatorname{Var}(x_i)=\lambda_i m_{2,A}+\sigma_{\epsilon,i}^2.
\]

After the flag-1 linear correction $y_i=a_i x_i$, the exact marginal
variance is

\[
\boxed{
\operatorname{Var}(y_i)
=a_i\kappa\operatorname{E}(y_i)
+a_i^2\sigma_{\epsilon,i}^2.
}
\]

The first term is compound-Poisson noise: it contains photon statistics and
the width of the MCP/phosphor/CCD pulse-height distribution. The second term
is additive electronic noise.

## Equivalent-output-quanta constraint

The instrument literature reports noise using equivalent output quanta,

\[
Q=\left(\frac{\operatorname{E}(Y)}
{\sqrt{\operatorname{Var}(Y)}}\right)^2.
\]

Therefore,

\[
\operatorname{Var}(Y)=\frac{\operatorname{E}(Y)^2}{Q}.
\]

For WIC, Mende et al. report $Q=23$ per 100 R exposure per resolution cell.
The resolution cell is $0.18^\circ$, while a WIC detector pixel is
$0.09^\circ$; one resolution cell therefore contains $2\times2=4$
detector pixels. Frey et al. report a reference sensitivity of 612.6 AD
units/kR per detector pixel during a ten-second exposure. At 100 R,

\[
\bar y_{\mathrm{pixel}}=0.1\times612.6=61.26\ \mathrm{AD},
\]

and the expected signal summed over one uniformly illuminated resolution cell
is

\[
\bar y_{\mathrm{cell}}=4\times61.26=245.04\ \mathrm{AD}.
\]

If additive noise is negligible at this level, the effective linear variance
coefficient at the reference calibration is

\[
\kappa_{\mathrm{eff},0}
=\frac{\bar y_{\mathrm{cell}}}{Q}
=\frac{245.04}{23}
=10.65\ \mathrm{AD}.
\]

The units provide a useful check: multiplying this AD-valued coefficient by
the saved signal in AD produces a variance in AD$^2$.

The corresponding reference-space candidate is

\[
\boxed{
\widehat{\operatorname{Var}}(y_i)
\simeq10.65\,y_i.
}
\]

This is $10.65$, not $61.26/23=2.66$, because the published 23 quanta
refer to a four-pixel resolution cell.

## Applying flat-field and mission corrections

The exact corrected-space signal-dependent coefficient is $a_i\kappa_i$.
The literature constrains this product only at the reference calibration. For
another pixel and date,

\[
\operatorname{Var}(y_i)
=a_i\kappa_i\operatorname{E}(y_i)
+a_i^2\sigma_{\epsilon,i}^2.
\]

If the raw pulse-height coefficient is assumed constant,

\[
\boxed{
\widehat{\operatorname{Var}}(y_i)
\simeq10.65\,a_i y_i.
}
\]

This is a correction-aware, conservative candidate rather than a uniquely
identified equation. A stellar sensitivity change can arise from at least two
mechanisms:

- If it is mainly a change in post-detection pulse gain, the mission correction
  also restores the pulse-height scale and much of its factor cancels from
  $a_i\kappa_i$. The reference law $10.65y_i$ is then closer.
- If it is mainly a change in photon-detection efficiency, the corrected image
  represents fewer actual detected events and the remaining factor $a_i$ is
  required.

The stellar correction table measures total sensitivity but does not separate
these mechanisms. A sensible validation should therefore treat
$10.65y_i$ and $10.65a_i y_i$ as lower and correction-aware candidates,
rather than conceal the distinction.

At zero and very low WIC signal, the direct plug-in estimate is not usable for
inverse-variance weighting. Treating the equivalent number of output quanta
with the Jeffreys-Poisson estimate adds half an equivalent event:

\[
v_i = \kappa y_i + \frac{1}{2}\kappa^2,
\]

where \(\kappa=10.65\) for the reference model. The resulting variance floor
is about 56.7 AD\(^2\), or a standard deviation of about 7.5 AD. For the
correction-aware candidate, replace \(\kappa\) by \(10.65a_i\). This prevents
zero-count pixels from being treated as perfectly known, but it does not
replace the still-unknown additive CCD/read-noise term.

## What the WIC approximation cuts

1. **Additive CCD noise.** The term $a_i^2\sigma_{\epsilon,i}^2$ is omitted.
   The available single-point $Q$ calibration cannot separately identify a
   signal-dependent slope and a constant intercept. This matters most at low
   signal.
2. **Sub-resolution-cell covariance.** The published $Q$ constrains the
   variance of a $2\times2$-pixel sum:

   \[
   \operatorname{Var}\!\left(\sum_{i=1}^{4}y_i\right)
   =\sum_i\operatorname{Var}(y_i)
   +2\sum_{i<j}\operatorname{Cov}(y_i,y_j).
   \]

   Assigning $10.65y_i$ independently to individual pixels assumes the
   covariance terms are zero and that a uniform signal is shared equally.
   The MCP spot, optical PSF, CCD charge footprint, and TDI mapping can violate
   that assumption.
3. **Uncertainty in the flat-field and mission correction.** These are treated
   as known multipliers. Their uncertainty would produce correlated,
   multiplicative errors rather than ordinary detector noise.
4. **Absolute sensitivity uncertainty.** The uncertainty in 612.6 AD/kR and
   the compatibility of the pre-flight $Q=23$ measurement with the in-flight
   day-171 stellar calibration are not propagated.
5. **Clipping and saturation.** The FUVVIEW3 expression clips the calibrated
   image at zero. A clipped value is not an unbiased compound-Poisson sample.
   Invalid and saturated pixels should be masked.
6. **Non-Gaussian low-signal behavior.** A Gaussian perturbation with the same
   variance is only an approximation to compound-Poisson statistics.

These limitations mean that the WIC equation is appropriate as a transparent
first-order input model and sensitivity test, not yet as a complete detector
covariance claim.

# The fuvpy WIC reflattening

By default, `fuvpy.read_idl()` applies `_reflatWIC()`. In simplified notation,
the operation is

\[
y_i'=\left(\frac{y_i}{\rho_i}-C\right)\rho_i+C
=y_i+C(1-\rho_i),
\]

where $\rho_i$ is fuvpy's two-segment WIC reflattening profile and $C$ is
estimated as a median over selected background pixels. The symbol $\rho_i$ is
used to distinguish this operation from the date-selected FUVVIEW3 calibration
factor $f_i$. If $C$
is held fixed, the derivative $\partial y_i'/\partial y_i=1$, so the marginal
detector variance is unchanged. Because $C$ is actually estimated from the
same image sequence, it also introduces shared covariance between pixels and
frames. The present first-order input model freezes $C$ and cuts this shared
term. It must be listed separately from the subsequent BS/SH background-fit
uncertainty.

# Propagation through a fitted background

Let the flattened image be $\mathbf y$, the background design matrix be
$\mathbf G$, the converged fuvpy row multipliers be
$\mathbf D=\operatorname{diag}(d_i)$, and the regularization matrix be
$\Lambda$. Conditional on the final mask, basis, weights, and
regularization, fuvpy solves

\[
\widehat{\mathbf m}
=\mathbf A^{-1}\mathbf G^T\mathbf D^2\mathbf y,
\qquad
\mathbf A=\mathbf G^T\mathbf D^2\mathbf G+\Lambda.
\]

For an input covariance $\mathbf R=\operatorname{Cov}(\mathbf y)$, the
sampling covariance of the fitted coefficients is

\[
\boxed{
\operatorname{Cov}(\widehat{\mathbf m})
=\mathbf A^{-1}\mathbf G^T\mathbf D^2
\mathbf R\mathbf D^2\mathbf G\mathbf A^{-T}.
}
\]

With the diagonal approximation
$\mathbf R=\operatorname{diag}(v_i)$, the middle matrix is

\[
\mathbf G^T\operatorname{diag}(d_i^4v_i)\mathbf G.
\]

The fitted-background covariance is

\[
\operatorname{Cov}(\widehat{\mathbf b})
=\mathbf G\operatorname{Cov}(\widehat{\mathbf m})\mathbf G^T.
\]

The corrected image is $\mathbf c=\mathbf y-\widehat{\mathbf b}$. Its
covariance is not the sum of the image and background variances because the
background was fitted from the same image:

\[
\boxed{
\operatorname{Cov}(\mathbf c)
=(\mathbf I-\mathbf H)\mathbf R(\mathbf I-\mathbf H)^T,
}
\]

where

\[
\mathbf H=\mathbf G\mathbf A^{-1}\mathbf G^T\mathbf D^2.
\]

For one diagonal element,

\[
\operatorname{Var}(c_i)
=v_i+\operatorname{Var}(\widehat b_i)
-2\operatorname{Cov}(y_i,\widehat b_i).
\]

The covariance term must be retained. Adding `input variance + fitted-model
variance` alone double-counts part of the input fluctuation.

This analytic expression freezes the converged robust weights. Complete
repeated fits are required to include changes in Tukey weights and, for WIC,
the dependence between the sequential BS and SH fits. The BS and SH component
variances must not be added as if the two stages were independent.

# Recommended working models

| Sensor | Input variance in saved-image units | Status |
|---|---|---|
| SI12 | $v_i=a_i y_i$, with an explicit zero-count rule | Opt-in loader integration implemented; complete-sequence validation remains |
| SI13 | $v_i=a_i y_i$, with an explicit zero-count rule | Opt-in loader integration implemented; direct integer-recovery verification is still desirable |
| WIC, reference | $v_i=10.65y_i+p(10.65)^2$ | Default opt-in WIC loader model; use $p=0.5$ for low-count weighting |
| WIC, correction-aware | $v_i=10.65a_i y_i+p(10.65a_i)^2$ | Implemented as `read_idl(..., wic_variance_model="correction_aware")`; sensitivity interpretation remains unresolved |

The current residual-derived laws, including
$v_{\mathrm{WIC}}=19.0424B$, must not replace these detector equations.
Those laws were obtained from background-fit residuals and include auroral
structure, spatial correlation, and model mismatch. They failed held-out
transferability for WIC and SI13. Likewise, `dgsigma` is a robust residual
scale used by the background algorithm; it is not the detector measurement
standard deviation.

# Terms deliberately outside this component

The proposed variances are conditional on the supplied calibration and the
chosen fuvpy background model. They exclude:

- uncertainty or bias in flat-field databases and mission sensitivity curves;
- absolute radiometric sensitivity;
- FUVVIEW3 or fuvpy background-model discrepancy;
- atmospheric dayglow and auroral-emission modelling;
- viewing geometry and line-of-sight corrections;
- WIC/SI coregistration, footprint allocation, and target-grid interpolation;
- proton correction and assumed proton energy;
- the WIC/SI13 electron-energy retrieval or Zhang--Paxton replacement;
- Robinson or other conductance forward-model uncertainty; and
- VAE model and sampling uncertainty.

These components may be propagated later, but combining them should occur in
the IMAGE-processing or Monte-Carlo pipeline, not by inflating the detector
variance until residuals happen to agree.

# References and evidence

The three papers below are stored in
`icBuilder/literature/`.

1. Mende, S. B., et al. (2000), *Far ultraviolet imaging from the IMAGE
   spacecraft. 1. System design*, Space Science Reviews, 91, 243--270.
   File: `Mende_et_al_2000_IMAGE_FUV_system_design.pdf`.
2. Mende, S. B., et al. (2000), IMAGE Spectrographic Imager instrument paper.
   File: `Mende_et_al_2000_IMAGE_SI.pdf`.
3. Frey, H. U., et al. (2003), *Summary of quantitative interpretation of
   IMAGE far ultraviolet auroral data*, Space Science Reviews, 109, 255--283.
   File: `Frey_et_al_2003_IMAGE_FUV_quantitative_interpretation.pdf`.
4. Recovered FUVVIEW3 source and support tables:
   `/home/bing/Downloads/fuview3/fuview3/`, particularly `getudf_var.pro`,
   `get_flatfield.pro`, `fuv_mission_correction.pro`, and `support/`.
5. Current fuvpy uncertainty implementation:
   `fuvpy/src/background.py`, `fuvpy/src/background_uncertainty.py`, and
   `fuvpy/src/load.py`.
