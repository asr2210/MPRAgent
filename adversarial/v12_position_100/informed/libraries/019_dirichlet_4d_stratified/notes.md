# 019 — 4D-stratified Dirichlet(0.3)

## What I built
500k Dirichlet(0.3) compositions, binned into 5×5×5×5=625 4D quantile bins, sampled
80 per bin. 313 bins were empty (4D corners that Dirichlet rarely visits).

## Result
- eval_01 = 0.0773 (vs Dirichlet(0.3) 0.0786). -0.0013 loss.
- mean ≈ 0.0948 (vs 0.0954). Essentially neutral.

## Interpretation
4D stratification did not break the ceiling. Pattern across stratification attempts:
- Random Dirichlet(0.3): 0.0786 (best)
- GC-stratified (014): 0.0775
- Mid-GC rejection (017): 0.0773
- 4D stratified (019): 0.0773

**All stratification hurts equally.** Natural Dirichlet(0.3) draws are nearly optimal.
Any deviation toward forced uniformity or restriction costs ~0.001-0.002 eval_01.

## Hypothesis killed
H19: 4D composition stratification breaks ceiling → REJECTED.

## What to try next
Test if MULTINOMIAL noise within sequence generation matters. Force exact composition
match (place exactly round(200·p) of each base, then shuffle positions).
