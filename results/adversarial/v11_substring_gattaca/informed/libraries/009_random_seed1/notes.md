# E9: random_uniform with seed=1 (variance check)

## Design
Identical to E2 but seed=1.

## Result
- eval_01: 0.8507 (vs E2 seed=0: 0.8565, **delta -0.006**)
- Mean across 14: 0.835 (vs E2 0.841, delta -0.006)
- SK-N-SH stable

## Interpretation — important noise calibration
Two seeds of random_uniform:
- seed=0 (E2): mean = 0.8408
- seed=1 (E9): mean = 0.8354
- Avg: 0.838, range: 0.006

**Single-seed run variance is ~0.006.** Effects smaller than 0.01 are
within noise; for a single-seed comparison I need a delta of at least
~0.01 to claim a meaningful improvement.

Notably E8 antithetic_rc landed at 0.839, also within noise of 0.838 baseline.

## Updated baselines
- random_uniform single-seed: 0.835-0.857 eval_01, mean 0.835-0.841
- Pooled estimate: 0.838 ± 0.005 SD

## Next
Anything that doesn't shift mean by ≥0.01 is meaningless under single-seed
budget. Need to either:
- Run multi-seed averaging (eats experiment budget)
- Try interventions that should produce larger effects (>0.02)

E10 candidate: contrastive mutation pairs — 25k random + 25k single-point
mutations of each. Per-position distribution preserved. If model benefits
from contrastive pairs (similar sequence, different labels), this could
produce a measurable lift.
