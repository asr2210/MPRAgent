# Experiment 002 — dhs_uniform_random

## Design
Pure uniform-random sample of 50k Meuleman 2020 DHS elements, 200 bp
centered on summit. Calibration check.

## Hypothesis
Should land near the published `dhs_random` = 0.7089 if my data pipeline
matches the baselines'. Should beat exp 001 (0.5600) if specificity
weighting was the bug.

## Result
`eval_01 = 0.5627` — marginally better than 001 (0.5600). But still ~0.15
below the instructions.md `dhs_random` baseline.

## Interpretation
Two reads:

1. **My weighting in 001 wasn't the main issue.** Uniform is only 0.003
   higher than the inverse-numsamples-weighted version. Both are far
   below the published 0.7089.

2. **The instructions.md baselines and the v13 evaluation are not on the
   same axis.** Looking at strategies.md (v13-specific), the listed
   baselines are *all synthetic*: `random_uniform = 0.1414`, `gc_sweep =
   0.4169`, `gc_50 = 0.1108`. None of the DHS / SEI / MPRA biological
   baselines from instructions.md appear. The best v13 synthetic
   baseline tops out at 0.42 — exactly the ballpark where uniform DHS
   should beat them by ~0.15, which matches what I see. **v13 looks like
   a stricter / different evaluation harness than the one that produced
   instructions.md's Table 1.**

So my 0.5627 is already ~0.14 above the best documented v13 baseline
(gc_sweep). I can stop chasing dhs_topic = 0.7232 and start asking what
pushes me above 0.57 within v13.

## Updated reference frame
- v13 floor (synthetic random): ~0.14
- v13 best synthetic baseline (gc_sweep): 0.42
- uniform DHS (my exp 002): **0.5627** (new biological floor)
- target: as high as possible above 0.5627

## Numbers
mean_r averaged across 14 evals: 0.516
eval_01: 0.5627 (cf. 0.5600 in exp 001)
eval_08: 0.1663 (still very low — synthetic eval, no synth in library)
time_s: 13.8

## Pattern across evals
All evals fall in 0.51–0.60 except eval_08 (0.17). The model trained on
genomic sequences extrapolates poorly to whatever eval_08 contains
(suspected to be synthetic / OOD). Mixing some non-genomic content into
the library should help eval_08 at a small cost elsewhere.
