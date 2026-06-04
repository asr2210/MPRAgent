# 015_pure_uniform_random — notes

## Design
50k purely uniform-random 200bp ACGT sequences. No genome.

## Result
eval_01 = 0.0420
**eval_08 = 0.1237  (BREAKTHROUGH — highest single-eval score seen)**
mean_r = 0.0453  (TIED with dhs_specific 0.0452!)

## Interpretation
KEY DISCOVERY: pure uniform random ties dhs_specific on mean_r because
eval_08 strongly rewards random/synthetic content (0.124 vs 0.04-0.06 for
bio libraries). All other evals modestly favor bio content.

The mean_r metric averages over evals with very different preferences.
A library that maximizes mean_r should COMBINE bio sequences (win on
13 evals) with synthetic random (win on eval_08).

## Next probes
- Exp 016: 50/50 dhs_specific + uniform — expected eval_01 ~0.046, eval_08 ~0.09
- Tune mixture ratio for max mean_r
- Try 33/33/33 dhs/uniform/cCRE
