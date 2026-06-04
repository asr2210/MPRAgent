# 009_promoter_tss

## Design
50k 200bp windows centered on RefSeq transcript TSSes (40,281 unique
TSSes on chr1–22+X+Y). Random ±50bp offset around each TSS. seed=0.
Motivated by Agarwal Nature 2024: 200bp promoter cores function as
cell-type-AGNOSTIC "on switches", potentially generalizing better
than enhancers or random tiles.

## Result — strong negative
- eval_01 = **0.3700** (vs random genome 002: 0.4992; Δ -0.129)
- eval_07 = 0.3849 (vs 0.5985; Δ -0.214) — huge drop
- eval_13 = 0.4188 (vs 0.6020; Δ -0.183) — huge drop
- eval_08 = 0.0816 (vs 0.0916; Δ -0.010)
- mean over 14 ≈ 0.379 (vs 0.503; Δ -0.124)

Promoter-only sequences are MUCH worse than broad random genome.
The "universal promoter" hypothesis fails in this pipeline.

## Interpretation
The published literature claims promoter cores are universal across
cell types. In this pipeline, training only on promoters cripples
generalization — model overfits to GC-rich CpG-island-like features
that don't match the broader distribution the eval is testing.

This generalizes the DHS finding (E3): any narrowly-curated
biological subset underperforms broad random genome. The narrower
the curation, the worse:
- E2 chr1/17/19/22 random: 0.4992 (broadest, best)
- E6 all-chrom random: 0.4823
- E3 DHS-stratified: 0.4378
- E9 promoter TSS: 0.3700 (narrowest, worst)

## Theory update v3

**Diversity beats focus in this pipeline.** The eval set's
distribution is broader than any single biological annotation
captures (promoter, enhancer, DHS). The model needs broad exposure
to learn the sequence prior the eval is testing.

Sub-hypothesis: chr1/17/19/22 (gene-rich) slightly beats all-chrom
because the eval has a mild gene-enrichment bias, but not full
promoter-restriction.

If true: the optimal library would be "random hg38 with light
gene-enrichment" — closer to E2 than to E9, E3.

## Implications for design

- Narrowly curated subsets (promoters, enhancers, DHS) are worse,
  not better.
- The eval set isn't testing how well a model fits any single
  functional element class; it's testing broad sequence prior.
- Subsequent experiments should EXPAND diversity in a meaningful way,
  not narrow it. Promising directions:
  - Cross-species real DNA (mouse): tests universal vertebrate prior.
  - Mixtures of random genome + cCRE + repeats: explicit element-class
    diversity within real DNA.
  - Activity-stratified across diverse sources using a fast predictor.

## What to try next: uniform DHS from full index
A quick test: pure dhs_random (uniform sample from full 3.6M DHS,
no stratification). In published baselines this beats dhs_stratified
by a small margin. In my pipeline this should be intermediate
between E3 (0.44) and E2 (0.50) if the "narrow biological focus"
hypothesis is right.

If E10 ≈ E2: DHS without stratification ≈ random genome — there's
no penalty for DHS-only IF you don't enforce stratification.
If E10 ≈ E3 (with stratification penalty being small/absent): DHS-
ness per se is what hurts.
