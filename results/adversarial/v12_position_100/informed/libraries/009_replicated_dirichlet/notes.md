# 009 — Replicated Dirichlet (noise averaging test)

## What I built
5000 unique Dirichlet(0.3) sequences × 10 copies each = 50,000 lines (shuffled).
Test: does the harness average duplicate (x,y) measurements? If yes, replication reduces
per-sequence noise by ~sqrt(10) ≈ 3×.

## Result
- eval_01 = 0.0732 (vs 0.0786 for exp 004 pure Dirichlet(0.3)) — **WORSE by 0.005**
- eval_08 = 0.0511 (vs 0.0716) — significantly worse on the highest-noise eval
- mean_r ≈ 0.0894 (vs 0.0954 for exp 004) — worse across the board

## Interpretation
**Replication did NOT help.** Two non-exclusive explanations:
1. The harness dedupes (only unique x → y pairs are used, so effective N = 5000 not 50k).
2. Per-sequence noise is not the bottleneck — sequence COVERAGE is. Losing 90% of unique
   sequences hurts model generalization more than noise-averaging helps.

The fact that eval_08 dropped the most (-29%) is suggestive: eval_08 is the noisiest eval
and would have benefited most from noise reduction if averaging worked. The opposite
happened. **Strong evidence the harness dedupes** (or weights uniques, not lines).

## Hypothesis killed
H9: noise averaging via replication is a lever → **REJECTED**

## What this tells me
- Future experiments should NOT replicate sequences.
- All 50k slots should be unique (or maximally distinct).
- The bottleneck at ~0.078 is NOT MPRA-level measurement noise.
- The bottleneck IS something about either (a) the composition variance ceiling within
  Dirichlet, (b) model-side underfitting that more diversity could help with, or (c) the
  3-cell-types-vs-14-eval coverage gap.

## What to try next (exp 010)
Refine Dirichlet alpha sweet spot. Current ladder: alpha=0.1 → 0.0752, 0.3 → 0.0786,
1.0 (implicit baseline) → 0.0765 (exp 002). Try **Dirichlet(0.5)** — interpolates between
0.3 (peak) and 1.0 (baseline). If 0.5 > 0.3, the peak is higher up; if 0.5 < 0.3, peak
is around 0.3 and the alpha lever is saturated.
