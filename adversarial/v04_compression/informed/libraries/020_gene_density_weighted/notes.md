# 020_gene_density_weighted

## Design
50k random 200bp tiles, sampled across all 24 hg38 chromosomes
weighted by LOCAL gene density (transcripts per 1Mb bin + 0.5 prior).
seed=0.

Gene-density top bins span chr17, chr6, chr19, chr15, chr16, chr11,
chr7, chr19 — much broader than my E2 chr1/17/19/22 cut.

## Result — CEILING BROKEN
- eval_01 = **0.5008** (vs E2 0.4992; Δ +0.0016)
- mean over 14 ≈ 0.504

This is the FIRST design above 0.50 in eval_01.

Per-eval improvements:
- eval_06: 0.5392 (vs E2 0.5285; Δ +0.011) ← biggest gain
- eval_03: 0.5404 (vs E2 0.5364; Δ +0.004)
- eval_07: 0.6005 (vs E2 0.5985; Δ +0.002)
- eval_10: 0.5399 (vs E2 0.5413; Δ -0.001)
- eval_13: 0.5968 (vs E2 0.6020; Δ -0.005)

Mixed but net positive. The biggest gain is on eval_06/11 (the same
eval); eval_13 slightly down.

## Interpretation
**Continuous gene-density weighting across ALL chromosomes beats the
hard chr1/17/19/22 filter.** Possible mechanisms:
1. Captures gene-rich tracts from chr11/16/6/15 that I missed.
2. Smoother gene-density distribution matches eval better.
3. The +0.5 epsilon prior provides regularization (still includes
   intergenic regions in proportion).

This validates a hypothesis from E15: chr19+22 alone narrowed too
much; the right strategy is "go gene-rich AND broad-chrom".

## Caveat — needs multi-seed verification

Δ +0.002 is just above the E14 noise floor (3-seed σ ~0.001).
Single-seed could be lucky. **E21 will repeat with 3 seeds.**

If 3-seed confirms ≥0.5005: real improvement. Push the angle.
If 3-seed drops to ~0.499: noise, treat as tied with E2.

## Theory update v9

The "0.50 ceiling" is BREAKABLE. The natural ceiling isn't broad
random; it's broad random WEIGHTED FOR GENE-RICHNESS, which adds
a tiny systematic lift by aligning the training distribution with
the eval distribution's mild gene-richness bias.

Updated mental model:
- Eval distribution: broad genome with mild gene-richness preference.
- E2 (chr1/17/19/22): captures part of the gene-richness, misses
  gene-rich tracts on other chromosomes.
- E20: captures more of the natural gene-rich distribution.
- The gain is small (~0.002) because chr1/17/19/22 was already
  90%+ of the way there.

## Plan for E21
3-seed of E20 design to confirm.

## Plan for E22+ (if E21 confirms)
- E22: even stronger gene-density weighting (count^2 instead of
  count + EPS). Tests monotonicity.
- E23: combine gene-density with other axes (e.g., 80% gene-density
  weighted + 20% top-DHS).
