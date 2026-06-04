# Experiment 006 — Dirichlet(0.5) + JASPAR motif insertion

## What I did
For each of 50,000 dirichlet(0.5) backbone sequences (200bp), inserted
2 motif instances sampled from JASPAR CORE 2024 vertebrates (879 motifs).
Each instance sampled from the motif's PFM, placed at non-overlapping
random positions.

## Result
- eval_01 = **0.1364** (vs dirichlet alone 0.1395 — WORSE by 0.003)
- All cell types slightly lower than pure dirichlet

## Interpretation
Adding explicit TF motif content **hurts**. The model does not appear to
use motif/PWM information from these sequences. Motif insertion is
just "noise" that displaces the diverse iid composition signal.

Combined with the exp 003 (Markov) result, the conclusion is now strong:
**this weak-model evaluator learns from per-sequence base composition
(and possibly low-order k-mer counts), NOT from positional motif content
or higher-order sequence structure.**

## Theory update — strong
The library design problem reduces to: **maximize composition diversity
across sequences while keeping each sequence iid (high complexity)**.
Things that help: dirichlet sampling, wide range of (pA, pC, pG, pT).
Things that hurt: motif insertion, repeats, Markov chains, real DHS.

## Next
Test whether **better simplex coverage** (Sobol/QMC) beats random
dirichlet sampling. If yes: small wins by improving composition coverage.
If no: dirichlet is near-optimal and I should look at orthogonal axes.
