# Experiment 014 — 007 exact replicate with SEED=2

## Design
Identical to 007 in every way except SEED=1 → SEED=2:
- 500k GC-50 random ACGT
- Score with Malinois
- Top 50k by max(K562, HepG2, SKNSH)

Selected sequence stats are nearly identical to 007:
- Mean Malinois selected score: K562=2.60, HepG2=2.36, SKNSH=2.62
- GC mean=0.511 std=0.034

So the SELECTION process is highly reproducible — the only difference
is which 500k random sequences were drawn, then which 50k were selected.

## Result — measures the noise floor
- eval_01 = **0.3936** (vs 007 = 0.3969 → **DOWN 0.0033**)
- mean_r = **0.3827** (vs 007 = 0.3861 → **DOWN 0.0034**)
- Per-cell: K562=0.617, HepG2=0.430, SKNSH=0.133

## Why this finding is critical
The single-seed difference between two IDENTICAL recipes is 0.003 on
mean_r. This is the noise floor for evaluating libraries with this
pipeline. It directly invalidates the "ranking" of many prior experiments:

| Exp | mean_r | Δ vs 007 | Within ±0.003 noise? |
|-----|--------|----------|----------------------|
| 007 | 0.3861 | 0      | reference |
| 014 | 0.3827 | -0.0034| BORDERLINE (the noise floor itself!) |
| 013 | 0.3845 | -0.0016| YES (noise) |
| 011 | 0.3839 | -0.0022| YES (noise) |
| 006 | 0.3860 | -0.0001| YES (noise) |
| 008 | 0.3835 | -0.0026| YES (noise) |
| 010 | 0.3818 | -0.0043| MARGINAL |
| 005 | 0.3805 | -0.0056| LIKELY signal |
| 012 | 0.3805 | -0.0056| LIKELY signal |
| 009 | 0.3774 | -0.0087| SIGNAL (worse) |
| 003 | 0.3653 | -0.0208| BIG SIGNAL (worse) — PLS-only |

The "best library is 007" claim was largely lucky-seed.
The 0.397 plateau is real, but the rankings WITHIN it are noise.

## What is signal, given noise=0.003?
- GC composition is signal: PLS-only (mean GC 0.62) → mean_r 0.365 (Δ=0.020 BAD)
- Real cCREs (009) → mean_r 0.377 (Δ=0.009 BAD) — Malinois bias on cCREs is real
- Within the GC-50 cluster, ALL strategies give 0.380-0.388, a 0.008 spread that
  is barely larger than noise.

## Revised theory v8
**The 0.397 plateau is a fixed-point property of the (model, training-size,
evaluator) tuple given the GC-50 random subspace as input distribution.**

To break the plateau, an experiment must:
1. Produce a >0.005 improvement (signal > 1.5σ above noise)
2. Move OUTSIDE the GC-50 random subspace in some structured way
3. Add information the random sampler can't capture

Within the subspace, oracle selection appears NOT to be a reliable lift —
the 007 vs random gap (~+0.005) collapses to ~+0.001 vs 014.

## What to try next
Need INTERVENTIONS with predicted effect size > 0.005:
1. **Motif-tiled sequences** (015): fill 200bp with 8 contiguous JASPAR motifs
   — fundamentally structured, NOT in GC-50 random subspace
2. **GC variance test**: per-sequence GC varies in [0.30, 0.70], library mean
   still 0.50 — tests per-sequence vs library composition
3. **Multi-pipeline-seed average**: pool selections from 5 different oracle
   seeds — does this stabilize and beat single-seed?
4. **Real cCRE no-GC-filter**: lifts GC restriction, gives natural composition
   distribution — different subspace
5. **Worst-of**: deliberately bad library (e.g., GC=0.30 only) to confirm
   downward range and check for asymmetry around 0.50

The 0.397 ceiling may be evaluator-structural. Need to probe whether
it's the SUBSPACE that's saturated, or the (model, evaluator) tuple.
