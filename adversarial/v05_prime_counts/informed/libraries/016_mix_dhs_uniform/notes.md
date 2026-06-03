# 016_mix_dhs_uniform — notes

## Design
25k dhs_specific + 25k uniform random (50/50).

## Result
eval_01 = 0.0469
eval_08 = 0.1035
**mean_r = 0.0473 — NEW HIGH** (vs previous best dhs_specific 0.0452)

## Interpretation
Mixing breaks the mean_r ceiling. Bio half satisfies eval_01-07/09-14;
uniform half drives eval_08 from 0.055 → 0.103. Linear interp on eval_08
would predict only ~0.090; we got higher → potential positive synergy.

## Next probes
- Exp 017: 12.5k dhs + 37.5k uniform (skew toward eval_08 favorite)
- Exp 018: 50k uniform but with seeded GC content matching DHS
- Exp 019: 4-source mix
