# 007 — Per-sequence Markov chain with Dirichlet(0.3) rows

## What I built
For each of 50k sequences: 4×4 transition matrix with each row ~ Dirichlet(0.3); generate 200bp from a random start. Per-sequence dinucleotide profile varies.

## Result
- eval_01 = 0.0765 — WORSE than pure Dirichlet(0.3) (0.0786). Basically tied with Dirichlet(1.0) baseline (0.0768).
- mean across evals: 0.0944 (vs 0.0954 in exp 004).

## Interpretation
Markov dinucleotide structure does NOT help on top of base-level composition. Two possible reasons:
1. Model can't learn dinucleotide axes at 50k (only learns 4-base composition).
2. Markov chain partially homogenizes per-sequence composition (stationary distribution from transition matrix is less variable than direct Dirichlet draws).

Most likely (2). The chain converges to a stationary base distribution determined by the transition matrix, and this stationary distribution is less extreme than directly drawing compositions from Dirichlet(0.3).

## What to try next
Test whether DHS biology becomes useful when paired with EXPLICIT compositional diversity within the DHS pool. EXP 008: GC-stratified DHS — bin Meuleman DHS by GC content into 10 bins, sample equal per bin. Forces compositional spread inside the biology pool. If this beats DHS-stratified (0.0739) AND beats Dirichlet(0.3) (0.0786), biology+composition together is the winning combination.
