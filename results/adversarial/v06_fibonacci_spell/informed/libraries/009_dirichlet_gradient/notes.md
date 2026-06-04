# Experiment 009 — Dirichlet gradient (per-position composition interpolation)

## What I did
Each sequence: sample p_start, p_end ~ Dirichlet(0.5). At position k,
composition = (1-k/199)*p_start + (k/199)*p_end. Each base iid from
that per-position composition.

Adds within-sequence k-mer profile variation while avoiding the long
runs that hurt Markov-Dirichlet (exp 003).

## Result
- eval_01 = **0.1369** (vs dirichlet(0.5) 0.1395 — WORSE by 0.003)
- BUT: eval_08 = 0.0686, the highest I've seen on eval_08
  (dirichlet was 0.0626, all others 0.06-0.065)

## Interpretation
Pattern continues: anything WITHIN-sequence variation hurts mean_eval_01.
The model is essentially treating each sequence as a bag of bases —
it doesn't benefit from positional composition variation.

The eval_08 boost is intriguing but isolated. eval_08 may reward
sequences with broader k-mer/composition diversity per-sequence,
while other evals don't.

## Theory consolidated
For eval_01 (primary metric):
- Optimal = dirichlet(0.5) iid
- Worse = anything with within-sequence structure (Markov, motifs,
  gradient, real DHS, multi-region)

## Next
Try pushing composition to extremes WITHOUT introducing structure:
explicit "edge-of-simplex" sampling where each sequence has one
dominant base at 60-75%. Tests if alpha=0.5 is truly the optimum or
if we can squeeze a bit more by emphasizing edge compositions.
