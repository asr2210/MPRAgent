# 005_satmut

## Design
100 seed cCRE sequences × 500 variants each (1-3 single-base substitutions).
Total: 50,000 sequences, but only 100 unique contexts.

## Hypothesis
If eval measures position-specific effects (per-base mutation impact), this
should score very differently. If eval is about overall regulatory grammar,
this should score modestly (only 100 contexts learned).

## Result
**eval_01 = 0.0046** — virtually zero. Many evals NEGATIVE.

## Interpretation
Diversity matters MORE than I appreciated. Going from 50k unique
genomic sequences (cCRE: 0.047) to 100 unique seed contexts × 500 variants
(satmut: 0.005) drops score 10×. And k562_r is NEGATIVE for all evals →
the model learned patterns ANTI-correlated with what eval measures.

This is a strong constraint: the library must have HIGH unique sequence
diversity. Any strategy that collapses to a small set of templates fails.

## What to try next
- Stop reducing diversity. Keep N_unique ≈ 50k.
- Try maximally-diverse libraries: multi-source mixes, random genomic
  windows from anywhere in hg38, etc.
- The path forward is likely about MATCHING the specific compositional /
  structural distribution of the eval, not about reducing variance.
