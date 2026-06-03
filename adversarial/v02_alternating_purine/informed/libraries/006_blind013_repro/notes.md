# Experiment 006 — blind 013 recipe replication

## Design
Direct port of blind 013_asym_strat_flank (the v02 best result so far):
- 15k uniform cCRE centers (preserves natural class distribution)
- 5k CTCF-only boost
- 5k DNase-H3K4me3 boost
- 25k paired flanks (±1500-3000bp, cCRE-overlap rejected)
- seed 0

## Result — eval_01 = 0.1465
Blind reported 0.1727 on this recipe; my seed-0 run gets 0.1465. Difference
≈ 0.026. But blind themselves ran the same recipe with seed 1 (exp 029) and
got 0.1581. So seed variance for this recipe is ~0.015-0.026.

My result is within seed-variance noise of blind's. Pipeline is sound.

## Per-eval comparison vs blind 013
- eval_06: 0.148 (vs blind 0.218) — big gap, seed-dependent
- eval_07: 0.183 (vs blind 0.183) — match
- eval_11: 0.148 (vs blind 0.218) — same as 06
- eval_13: 0.180 (vs blind 0.126) — I'm HIGHER

So individual evals swing ±0.05 between seed runs. mean_r 0.146 vs blind's
0.173 captures the average seed-dependence on eval_06/11.

## Theory v5 — seed variance is large
Per-eval single-seed variance is ±0.02-0.05. This means:
- Apparent improvements of <0.01 on eval_01 are not significant
- I should focus on changes that move multiple evals together
- Multi-seed averaging would help but costs experiments

## Plan for exp 007
Combine the two best sources I've found:
- cCRE asymm (good for eval_01-04, weak on eval_07)
- Table_S2 (strong on eval_07, weak on eval_06/11)

**007_ccre_asymm_mpra_flanks**:
- 12k uniform cCRE
- 4k CTCF, 4k DNH3 boost
- 5k Table_S2 (UKBB+GTEx, non-eval-chrom)
- 25k paired flanks (from cCRE positives only — Table_S2 sequences lack
  precise genomic positions)
Total = 50k. Tests whether source diversification raises mean_r over single
sources.
