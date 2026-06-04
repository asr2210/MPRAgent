# E014 — Gosai chrom-balanced (23 chroms, ~2174 each)

50K Gosai sequences with equal representation per chromosome (chr 1-22 + X).

## Result
eval_01 = 0.3218. Per-cell K562=0.149, HepG2=0.203, SKNSH=0.613.

## Interpretation
Chrom-balance ≈ random Gosai (0.323) — slightly worse. This rules out
chrom-balance as the lever behind E008's boost. E008 (chr 7,9,13,21,X
only) at 0.336 is genuinely BETTER than chrom-balanced 0.322.

**Strong evidence for eval-leak**: those 5 specific chroms appear to
overlap eval-set sequences. Time to confirm by running the OPPOSITE:
Gosai excluding those chroms.
