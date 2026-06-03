# 008_genome_aug_motif

## Design
50k real-genome tiles (chr1/17/19/22, same as 002) with one strong
TF consensus motif planted at a random position per sequence.
Motif library: 53 entries, random strand. seed=0.

## Result
- eval_01 = **0.4958** (vs 002's 0.4992; Δ -0.003)
- mean over 14 ≈ 0.495
- All evals essentially identical to 002 within ±0.005.

## Interpretation — null result
Augmenting real DNA with ectopic strong motifs does NOT improve over
plain real DNA. Real DNA is not motif-content-limited; the model
already extracts what it can from natural motif occurrences in
context. Ectopic motifs in non-natural positions add no signal.

## Implications
- The ceiling for 50k libraries in this pipeline appears to be ~0.50
  for eval_01.
- Brute-force motif augmentation is not the path forward. Need a
  qualitatively different design axis.
- Possible directions for the ceiling:
  - Diversity expansion via cross-species / curated-element sources
  - Activity-stratified design (needs an oracle predictor)
  - Mixed-element design: many different functional element classes
    represented (cCRE: promoter / enhancer / CTCF / etc.)
  - Library skewed to the eval distribution itself (if I could
    characterize it)

## Theory update v2
After 8 experiments, the picture:
- Single biggest effect: real DNA vs synthetic random (+0.19).
- Within real DNA, refinements via DHS-enrichment, broader chromosomes,
  motif augmentation are all small or negative (-0.05 to +0.00).
- 4-mer composition matching is *anti-productive* — sequences that look
  like real DNA without structure confuse the model worse than featureless
  random.
- Isolated motif content in random scaffolds recovers only ~20% of the
  real-DNA prior.

The real-DNA prior is dominated by something I haven't isolated:
likely the combination of repeat structure, gene-body / intergenic
context, and natural motif co-occurrence grammar. This is hard to
synthesize and hard to engineer; real DNA is currently the only
known source.

Next: read the literature on what's been tried for MPRA library
design before another experiment.
