# Experiment 003 — Markov-Dirichlet (first-order Markov chains, per-seq Dirichlet)

## What I did
Each sequence = first-order Markov chain with its own 4x4 transition matrix
(each row ~ Dirichlet(0.5)) and own initial distribution. Adds dinucleotide
diversity on top of base composition diversity.

## Result
- eval_01 = **0.1346** (worse than 002 dirichlet 0.1395 by 0.005)
- Better than 001 synthseqs (0.1319) by only 0.003
- K562 = 0.035 (between synthseqs 0.029 and dirichlet 0.047)

## Why this hurt
Markov chains can produce low-complexity sequences (e.g., homopolymers
if the transition matrix has a high diagonal). The v06 baseline shows
homopolymer_rich = 0.058 (very bad) — low-complexity sequences are
actively harmful. By varying transition matrices via Dirichlet, I
inadvertently created a fraction of low-complexity sequences.

Inspection of generated seqs confirms: many begin with long runs like
"CCCCCCG..." or "GTTTGTTT...".

## Theory update
i.i.d. diverse compositions > Markov diverse compositions. The lesson:
**diversity in base composition between sequences helps, but introducing
structural/dinucleotide-level patterns WITHIN sequences hurts**, because
the model gets sequences that look "non-random" in ways that don't
correspond to real regulatory grammar.

Generalization rationale unchanged: i.i.d. sequences span composition
space cleanly without confusing the model with spurious dinucleotide
correlations.

## Next
Test whether **real sequences add value when mixed with the best
synthetic strategy** (dirichlet). If yes, biology matters. If no,
the diversity story dominates.
