# 004_dhs_genome_mix

## Design
25k random hg38 tiles (chr1/17/19/22, same source as 002) + 25k DHS-
stratified sequences (same source as 003). Total 50k. seed=0.

Direct test of DHS↔genome complementarity.

## Result
- eval_01 = **0.4926** (vs 002: 0.4992; 003: 0.4378)
- eval_07 = 0.5956 (vs 002: 0.5985; 003: 0.4521)
- eval_13 = 0.5877 (vs 002: 0.6020; 003: 0.4672)
- eval_08 = 0.0917 (vs 002: 0.0916; 003: 0.1016)
- mean over 14 ≈ 0.495 (vs 002: 0.503; 003: 0.439)

## What this means
Naive linear interpolation between 002 and 003 predicts 0.469. We got
0.495. So the mix is *better than the mean of its parts* — there is
some complementarity. But it's still ~0.007 below pure genome-random
(within seed-noise).

Operational read: at 50k slots, replacing 25k genome sequences with
25k DHS sequences costs ~0.007 r-points. DHS sequences are slightly
less informative per slot, but not catastrophic. A focused-mixture
design is unlikely to beat pure broad-genome by a meaningful margin
unless I find a genuinely complementary source.

## Theory update
The "biological enrichment" axis (DHS-vs-random within-genome) is
not where the win comes from in this pipeline. The big effect is
"real DNA vs synthetic" (Δ +0.19). Within real DNA, the additional
selection on accessibility/NMF programs is a small effect at best.

So my generalization theory needs to shift away from "expose the
model to functional grammars" and toward "expose the model to the
broad genomic substrate distribution". A model that learns the
generic DNA prior (CpG, repeats, gene bodies, intergenic, etc.) and
trains on enough varied sequences will generalize. Cell-type-specific
chromatin annotation is not a useful axis for selecting sequences.

## What to try next
Two directions:
1. **Confirm "broad genome is the recipe"** — random tiling from a
   broader chromosome set (not just chr1/17/19/22). If similar to E2,
   robust. If much higher, broaden further.
2. **Test if a different functional axis adds**: activity diversity
   (sequences predicted to span low/mid/high MPRA score), conservation
   (PhyloP-weighted), or motif-implanted scaffolds.

The first is the more decisive next step because it tells me whether
my "real DNA is the recipe" finding generalizes beyond my specific
chromosomes.
