# Handoff - Latest

Last updated: 2026-09-08
Branch: `background_model_revision`
Checkpoint: Step 1 complete

## Completed

The released `backgroundmodel_BS` now has an explicit
`noise_model="binned"` option. It retains the legacy dayglow coordinate and
all other released fit behavior, but corrects the empirical `dgsigma`
regression by fitting one point per populated geometry bin. The default
`noise_model="legacy"` preserves the released result for comparison during
combined validation.

The bin edges define one fitted nightside bin below zero and 14 dayside bins
from 0 to 3.5 in steps of 0.25. Within each populated bin the calculation
forms the weighted means of `dgmodel**2` and `(img-dgmodel)**2`. A bounded
linear fit constrains both slope and intercept to be nonnegative, and
`dgsigma` is the square root of the fitted variance at every pixel.

## Verification

`pytest -q tests/test_background_dgsigma.py` passes 3 tests. The tests confirm
the required bin edges, recover known slope and intercept behavior, and show
that duplicating observations within one bin neither adds regression points
nor changes the fitted spread.

## Next action

Implement Step 2 from `vault/02_Algorithm/Background Model Revision
Implementation Plan.md` in a separate commit. Do not add experimental
geometry, uncertainty, diagnostics, generated figures, or unrelated
abstractions.

## Portfolio impact

The clean background-model revision has completed its first scientific
change. This is an implementation checkpoint only; package defaults and
downstream production guidance are unchanged pending combined validation.
