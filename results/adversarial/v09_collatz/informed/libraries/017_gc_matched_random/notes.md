# E017 — GC-content matched random uniform

For each of 50K sequences, sample a GC fraction from Gosai's GC distribution
(mean 0.462, std 0.107), then generate a random sequence with that GC.

## Result
eval_01 = 0.3059. Per-cell K562=0.137, HepG2=0.190, SKNSH=0.590.

## Interpretation — striking
Essentially identical to E016 shuffled Gosai (0.3052). So:

**Per-sequence GC distribution alone explains ALL of the +0.06 lift
from uniform random to shuffled Gosai.**

This means the pipeline's eval correlation is dominated by the model
learning to associate GC content with activity. Higher-order structure
(A vs T asymmetry, dinucleotide patterns, motifs, grammar) collectively
contribute only ~+0.018.

## Implication for ceiling
The 0.34 ceiling decomposes as roughly:
- Baseline (uniform random)        ≈ 0.24
- + Match Gosai GC distribution    ≈ 0.30 (+0.06)
- + Real grammar/motifs            ≈ 0.32 (+0.02)
- + Test-chrom eval-leak           ≈ 0.336 (+0.013)
- Hard ceiling at                  ≈ 0.34
