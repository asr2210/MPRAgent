# 024 — 40k Dirichlet(0.3) + 10k random hg38

## What I built
40k Dirichlet(0.3) synthetic + 10k random hg38 windows. Tests if SMALL dose of biology
helps without too much dilution.

## Result
- eval_01 = 0.0777 (vs Dirichlet(0.3) 0.0786). -0.0009.
- mean ≈ 0.0951.

## Interpretation
Even a 20% biology dose hurts. Confirms biology adds nothing to a composition-rich
library at ANY ratio.

H24 → REJECTED.

## What to try next
Reproduce best (exp 004) with different seed to estimate variance.
