# 013_human_topsignal_dhs

## Design
40k random tiles from human chr1/17/19/22 + 10k DHS sites from top 1%
mean_signal across 733 Meuleman biosamples (signal ≥ 3.70). Shuffled.
seed=0.

## Result
- eval_01 = **0.4956** (vs E2 0.4992; Δ -0.004) — null
- mean over 14 ≈ 0.494
- All evals within ±0.005 of E2.

## Interpretation
20% replacement with "verified universal regulators" makes no
difference. The hypothesis that top-mean_signal DHS sites would
add cell-type-transferable signal beyond random genome is FALSE.

Three negative-mix results now:
- E4 (50/50 DHS-strat + random): 0.4926 (Δ -0.007)
- E12 (50/50 mouse + human): 0.4839 (Δ -0.015)
- E13 (20/80 top-DHS + random): 0.4956 (Δ -0.004)

Smaller fraction = smaller penalty. There's no diversity bonus to
balance the dilution. Adding ANY non-random-genome content to pure
human random hurts (slightly) or is neutral.

## Theory update v6
**Random tiles from gene-rich human chromosomes is a strong local
optimum.** I cannot beat 0.50 by:
- Adding curated DHS at any fraction (E3, E4, E10, E13)
- Adding mouse cross-species (E11, E12)
- Adding motifs to random (E7) or to genome (E8)
- Replacing with Markov-synthesized real-looking DNA (E5)
- Narrowing to functional annotations (E9 promoter)

What hasn't been tried:
1. **Quality filtering**: reject simple-repeat / low-entropy tiles
   from random; tests if median random has wasted slots.
2. **Wider gene-rich chromosome set**: chr1/11/16/17/19/22 (6 chroms);
   tests if "even more gene-richness" helps.
3. **Activity-stratified across the full mean_signal range**: equal
   slots per signal decile (combines diversity with broad coverage).
4. **Mix of element CLASSES from cCRE annotations** (promoter,
   distal enhancer, CTCF-bound, intergenic) — different cuts than DHS.
5. **Multi-seed** — am I being fooled by single-seed variance?

## What to try next
**E14: human + activity-stratified DHS (10 deciles)** — 35k random
human + 15k DHS sites with 1.5k from each mean_signal decile of the
DHS index. Combines broad-genome backbone with explicit activity-
diverse DHS coverage (not just high-signal as in E13).

If this also ≈ 0.49-0.50: pretty conclusive that no DHS injection
breaks the human-random ceiling at any fraction. I'd then pivot to
quality-filtered random (E15) and multi-seed variance check.
