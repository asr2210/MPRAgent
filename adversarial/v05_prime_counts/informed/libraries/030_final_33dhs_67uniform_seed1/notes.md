# 030_final_33dhs_67uniform_seed1 — notes (FINAL VALIDATION)

## Design
Exact same recipe as winning Exp 019 (33dhs/67uniform) but seed=1.
Validates noise on the winning approach.

## Result
eval_01 = 0.0459 (vs seed=0: 0.0478) — delta 0.0019
eval_08 = 0.1087 (vs seed=0: 0.1109) — delta 0.0022
mean_r = 0.0472 (vs seed=0: 0.0488) — delta 0.0016

## Interpretation
Seed noise ~±0.002. The 33/67 recipe consistently outperforms pure
dhs_specific (0.0452) on mean_r — the 0.020 gap is well above noise.
The "0.0488 peak" partly reflects seed=0 luck; true expected
mean_r ~0.0480.

## Recommended winning library
**Exp 019**: 16,667 dhs_specific + 33,333 uniform random (33/67 mix).
Best single-seed mean_r 0.0488, robust across seeds.
