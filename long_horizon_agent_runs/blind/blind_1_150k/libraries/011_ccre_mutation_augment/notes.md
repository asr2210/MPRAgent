# 011 — cCRE + mutated cCRE pairs (10% substitution)

## Goal
Test whether data augmentation via small base perturbations (10 % per
sequence) gives the model contrastive training signal that improves over
150k unique cCRE sequences alone.

## Method
- Sample 75,000 unique cCRE-centered windows (same pool & filter), seed=10
- For each, generate a mutated version: 20 positions (10 %) substituted to
  a different ACGT base (guaranteed change), seed=10
- Library = 75k originals + 75k mutated = 150k. Shuffle.

## Result: regression. Mean = 0.8837 (vs 003 pure cCRE 0.8862, **−0.003**).

| eval | cCRE | mut-aug | Δ |
|------|------|---------|---|
| 01 | 0.8288 | 0.8253 | −0.003 |
| 02 | 0.9270 | 0.9237 | −0.003 |
| 03 | 0.9193 | 0.9151 | −0.004 |
| 04 | 0.8657 | 0.8629 | −0.003 |
| 05 | 0.8285 | 0.8251 | −0.003 |
| 06 | 0.9274 | 0.9240 | −0.003 |
| 07 | 0.9025 | 0.9021 | 0.000 |
| 08 | 0.8922 | 0.8974 | **+0.005** |
| 09 | 0.9465 | 0.9426 | −0.004 |
| 10 | 0.9272 | 0.9255 | −0.002 |
| 11 | 0.8145 | 0.8113 | −0.003 |
| 12 | 0.7959 | 0.7924 | −0.003 |
| 13 | 0.9038 | 0.9008 | −0.003 |
| 14 | 0.9276 | 0.9242 | −0.003 |

## Key observations
1. **Mutation augmentation alone doesn't help.** Net −0.003 mean.
2. **Eval_08 gains +0.005** — mutated cCREs act like a *weak* synthetic
   diversity component (out-of-distribution enough to trigger a little
   eval_08 bonus, but not as strong as motif-embedded random which gave
   +0.017 in exp 006).
3. **All other evals lose ~−0.003.** This is the cCRE diversity cost
   (75k unique vs 150k unique in exp 003).
4. **Hard evals 11/12 unchanged** (well, regressed by 0.003 like
   everything else).

## Theory update (v10 → v11)
- Mutation augmentation creates pairs of (sequence, activity) and
  (mutated, activity'). The model could in principle learn from the
  paired examples — but at 10 % mutation rate, the perturbations are
  small and the activity changes are tiny, so each pair adds little
  unique information.
- The "expand the cCRE pool" tactic continues to lose:
    * 007 class-balanced: −0.004
    * 008 non-cCRE genome diversity (10 %): −0.003 vs pure cCRE
    * 010 neighbourhoods: −0.004
    * 011 mutation pairs: −0.003
  All four lose because they reduce *unique-cCRE* coverage in exchange
  for something the model doesn't extract value from.
- The mutated cCRE gives a small eval_08 bonus, consistent with
  "any non-cCRE-distribution sequence helps eval_08 a bit". But the
  bonus is small enough that it can't offset the cCRE-diversity cost.

## Next
Going with experiment 012 = **reverse-complement augmentation** combined
with the winning hybrid:
- 67,500 unique cCRE-centered windows (forward)
- 67,500 same cCREs in reverse complement
- 15,000 motif-embedded synthetic
- Total 150k

Tests whether RC augmentation gives the model strand-aware training
signal that lifts performance, *and* whether the synthetic eval_08 bonus
survives at slightly reduced unique-cCRE count. Direct comparison to
exp 006 (mean 0.888, 135k unique cCREs + 15k motif).

Predicted: if RC helps, ≥ 0.888; if not, drops back to ~0.884.
