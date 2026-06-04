# Experiment 012 — 4-way hybrid mix

## Design
50k sequences = 12.5k from each of:
1. Random + Malinois top (from 007 pool)
2. Random + Malinois span (from 006 pool)
3. cCRE GC-matched (fresh sample)
4. GC-50 i.i.d. random (fresh)

All GC ≈ 0.50; final mean GC 0.501, std 0.035.

## Hypothesis
Mixing complementary distributions might exceed any single strategy if
diverse training signal helps the model generalize.

## Result
- eval_01 = **0.3903** (vs 007 = 0.3969)
- Mean across 14 evals = **0.3805** (vs 007 = 0.3861)
- Per-cell-type on eval_01: K562 0.607, HepG2 0.428, SK-N-SH 0.136
- Runtime: 1189s

## Interpretation
**Mixed strategy underperforms best single strategy.** Dilution effect:
12.5k of the best (007) plus 37.5k of weaker libraries averages down to
the weighted-average of those libraries' standalone scores, roughly.

This says: the model TRAINS its weights based on the FULL training set,
and mixing weak training distributions with strong ones dilutes the
strong signal rather than complementing it. There is no "diversity bonus"
at this scale.

## What this rules out
- "Strategy diversity helps generalization" hypothesis
- Mixing different sub-distributions as a way to beat the ceiling
- The model's ability to extract complementary value from heterogeneous
  training data

## Theory v6 reinforced
The 0.397 ceiling holds. The best library is the most internally
homogeneous, oracle-selected one (007). Mixing hurts.

## What to try next
- BIGGER candidate pool for 007's strategy (4-10× larger random pool to
  find more extreme top activity)
- Replicate 007 with different seed to measure score-noise
- Pure-composition extremes: exactly 100/200 GC bases per sequence
- Multi-oracle ensemble (Malinois + train a small CNN on the 007 library)
