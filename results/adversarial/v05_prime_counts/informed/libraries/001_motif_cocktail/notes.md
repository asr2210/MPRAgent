# 001_motif_cocktail

## Design
Each of 50,000 200bp sequences = random ~50% GC background with 4-8 JASPAR 2024
CORE motifs (sampled from PFM probabilities, length 6-20 each) planted at random
non-overlapping positions. 2258 motifs available.

## Hypothesis
Planted motif content should beat pure-random synthetic (synth_oracle eval_01
= 0.684 in instructions.md baselines). If it beats DHS (dhs_topic = 0.723),
motif density > genomic context.

## Result
**eval_01 = 0.0375** — essentially zero.
All evals near zero; eval_08 highest at 0.111. Pattern matches strategies.md
`random_uniform` (0.0399) and `gc_50` (0.0397).

## Interpretation
Motif-cocktail performance ≈ random_uniform. Two possible takeaways:

1. **The instructions.md baseline table does not predict this prepare.py's
   behavior.** In that table synth_oracle gets 0.684; here random_uniform
   gets 0.04. The two are inconsistent. The strategies.md (zero-ish for
   synthetic) is the reliable signal for this environment.
2. **Planted-motif synthetic sequences offer no learnable signal in this
   environment.** The model needs real genomic context to learn anything.
   Motif densities in random backgrounds may even be unrealistically high
   (4-8 motifs / 200bp), or the surrounding random bases lack the
   correlational structure the eval data uses.

## What to try next
- **Exp 002**: Real DHS sequences sampled with topic weighting (matches the
  best baseline `dhs_topic`). This is the highest-confidence next move — we
  need genomic context, period.
- Need: hg38 reference (downloading in background) + DHS index (downloaded).
