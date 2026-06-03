# Experiment 011: Markov-1 random with human dinucleotide frequencies

## Result
- eval_01: **0.4873** (mean_r)
- K562: 0.9131, HepG2: 0.5533, SKNSH: -0.0046

## Comparison
- random_seed42 (best): 0.5221
- Markov-1 (this):       0.4873
- **DROP of -0.035 vs i.i.d. uniform random**

## Interpretation
Natural dinucleotide context (esp. CpG suppression at C→G = 0.022) HURTS the model.
- K562 r dropped from ~0.99 to 0.913 — adding natural-feeling local correlations
  triggers the K562 oracle's penalty, the same way real biological sequences do.
- HepG2 r barely changed (0.553 vs ~0.56) — confirms HepG2 is invariant to library.
- SKNSH r remains ~0 — confirms no SK-N-SH signal recoverable from random-like data.

## Conclusion
- Even subtle natural structure (dinucleotide frequencies) is read as "biological"
  by the K562 oracle and degrades the score.
- The optimum is genuinely *unstructured* random, not random-shaped-like-genome.
- This eliminates "natural local context" as a useful intervention.
