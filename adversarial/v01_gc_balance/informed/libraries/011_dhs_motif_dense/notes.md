# 011_dhs_motif_dense — notes

**Design**: PWM-scan a 200k DHS subsample from the non-label-aligned pool with 80 JASPAR motifs (strided from 879). For each DHS, count motifs with at least one hit ≥80% of the motif's max log-odds score. Sample top 50k by this "distinct strong motifs present" score.

Pool mean = 3.51 motif hits/seq. Selected top-50k mean = 7.15. Selected DHS range: [5, 17] hits.

**Result**: eval_01 = 0.6084 (vs exp 005 = 0.6752, **-0.067**). Mean across 14 evals = 0.5683 (vs 0.6282, -0.060). Drop is broad: eval_04 falls hardest (0.3637 vs 0.5374). eval_07/13 ~unchanged.

**Motif-density selection HURTS.** Strongly.

**Why**: Selecting motif-rich DHSs CONCENTRATES the library on a particular kind of sequence (likely TF-cluster regulatory hotspots, e.g., super-enhancer-like regions with many adjacent TFBSs). This loses the background-context diversity of "typical" regulatory regions. Same failure mode as exp 007 (numsamples weighting).

**Theory update H4 → H5**:
> **Concentration on ANY single quality/property axis hurts.** Whether you weight by biological robustness (numsamples), TF-binding density (motif scan), or class balance (forced topic balance), pushing the library toward a narrow subset reduces the natural variance the model relies on. The optimal strategy is *uniform sampling from a thoughtfully filtered pool* (e.g., exclude label-aligned topics), not active selection within that pool.

This is a STRONG, falsifiable theory: any single-axis concentration should hurt. Counter-test would be motif-SPARSE selection — also should hurt (predicts ≈0.60 not 0.68).

**Implication**: stop "selecting" DHSs. Start *augmenting* them. Per-DHS multi-view (RC, offset) is the orthogonal axis: same biological pool, more views per source.

**Next (012)**: reverse-complement augmentation. Take 25k DHSs from the no-label pool, include each as forward + RC = 50k. Doubles strand coverage from the same biology. If conv net already learns RC-equivariance, this is flat (like positional jitter was). If it doesn't, this is a free lift.
