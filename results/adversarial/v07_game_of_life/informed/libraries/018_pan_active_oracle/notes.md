# Experiment 018 — Pan-active oracle (top by MIN cross-cell)

## Design
500k GC-50 random, score with Malinois, top 50k by MIN(K562, HepG2, SKNSH).
Selects sequences active across ALL THREE training cell types.

## Selected stats
- min-cell score: mean=2.17 (selected min over 3 cells per seq)
- K562 mean=2.56, HepG2 mean=2.41, SKNSH mean=2.70 (all comparable to 007)
- The MIN-active selection gives sequences with high baseline activity in
  every cell — vs 007's MAX selection where sequences have high activity in
  ONE cell.

## Result
- eval_01 = **0.3971** (Δ vs 007 = +0.0002)
- mean_r = **0.3863** (Δ vs 007 = +0.0002)
- Per-cell: K562=0.618, HepG2=0.435, SK-N-SH=0.138

## Interpretation
**Pan-active criterion matches max-active criterion exactly.** Within noise
of 007/017.

**Theory v12.** The oracle SELECTION CRITERION (max, min, span, discrim) is
irrelevant — they all hit the same ceiling. What matters is:
- The sequences are GC-50 random (natural-like)
- They are NOT extreme outliers (014 vs 013 confirmed)
- Selected from the active tail (not bottom, not middle — see 022)

The (model + 50k + evaluator) tuple cares about training-data being
natural-like and informative, not about WHICH information is selected
from within the active set.

## What this rules out
- "Pan-active sequences transfer better" hypothesis
- "Criterion choice within active sequences" as a lever

## What to try next
The oracle-selection axis is now exhausted. Remaining axes:
1. Real cCREs unfiltered (natural composition variation) — 019
2. RC-augmented library — 020
3. Bottom-active anti-control — 022
4. K-mer matched to cCREs — diff structural sequence statistics
5. Iterative oracle / self-training
