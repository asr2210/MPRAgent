# Experiment 021 — Mutation-around-oracle (5k seeds × 10 mutations each)

## Design
Take first 5,000 sequences from 007's selections (oracle-active GC-50).
For each, generate 9 single-base-mutated variants (different random
position each time, different random base). Total = 5,000 × 10 = 50,000
unique sequences, but with strong intra-cluster correlation.

GC: 0.510 ± 0.033 (unchanged from 007), 50k unique sequences.

## Result
- eval_01 = **0.3706** (Δ vs 007 = -0.026, SIGNIFICANT)
- mean_r = **0.3600** (Δ vs 007 = -0.026, SIGNIFICANT)
- Per-cell: K562=0.597, HepG2=0.409, SK-N-SH=0.106 — all drop

This is one of the biggest negative deltas in our sweep.

## Critical reconciliation with 020
At first glance, this CONTRADICTS exp 020 (25k unique + 25k RC ≈ 007).
But the resolution is:
- 020's 25k unique are INDEPENDENT (each i.i.d. random)
- 021's 5k seeds × 10 mutants are HIGHLY CORRELATED (each cluster shares
  199/200 positions)
- The effective info content of 021 << that of 020 << 007

So 020 is NOT saying "any 25k will do" — it's saying "25k INDEPENDENT
samples + their RCs is enough." 021 has 5k effective independent samples,
which is insufficient.

## Theory v14
The bottleneck is INDEPENDENT samples from the natural-like subspace.
With:
- 5k independent + 9 mutants each: mean_r 0.360
- 25k independent + 25k RCs: mean_r 0.386 (≈ 50k independent)
- 50k independent random: mean_r 0.386 (007/017/020)

So the gradient is roughly:
| Independent samples | mean_r |
|--------------------|--------|
| 5k | 0.360 |
| 25k | 0.386 |
| 50k | 0.386 |

Diminishing returns set in around 25k independent samples. Below that,
each missing independent sample costs ~0.001 mean_r (extrapolating).

## What this rules out
- "5k seeds × 10 augmentation" as a viable info-density strategy
- Active-set neighborhood exploration as a way to compress info
- Mutation augmentation as a free lunch

## What this opens
1. **Test 10k seeds × 5 mutations** (would give ~10k effective independent)
   — predict mean_r ~0.375 (interpolating).
2. **Test 25k seeds × 2 mutations** — should give ~25k effective ≈ 0.386
3. **In other words: independent sample count is the lever; if mutating,
   one needs at least ~25k seeds.**

## What to try next
1. **Bottom oracle (anti-control)** (022): top 50k by MIN Malinois —
   directly tests oracle directionality
2. **Independent samples scaling** — but we've shown 50k saturates
3. **Multi-oracle ensemble** — different prediction targets
