# Experiment 004 — Dirichlet + SynthSeqs 50/50 mix

## What I did
25,000 Dirichlet(0.5) synthetic + 25,000 SynthSeqs (topic-weighted real DHS).
Shuffled and saved.

## Result
- eval_01 = **0.1359**
- Dirichlet alone (002): 0.1395
- SynthSeqs alone (001): 0.1319
- Mix landed BETWEEN, but below the better individual.

## Interpretation
Adding real DHS to diverse synthetic HURTS overall performance — the
combined library underperforms 50k dirichlet alone. Half-sized dirichlet
(25k seqs) gives less signal than full 50k, and the 25k DHS doesn't
compensate.

**Real DHS sequences are dead weight in this evaluator.** Under the v06
prepare.py model, biological realism contributes nothing beyond what
diverse synthetic provides.

## Theory update
Confirmed: diversity > biology. For the rest of the run, I should
focus on PUSHING SYNTHETIC DIVERSITY rather than mixing in real seqs.

Caveat for generalization: this is specific to this weak evaluator.
A stronger downstream model (CNN, transformer) trained on real sequences
might extract motif grammar that dirichlet sequences can't teach.
But for THIS evaluation, synthetic diversity dominates.

## Next
Push composition diversity further. Options:
- Smaller alpha (more extreme compositions per sequence)
- Mixture of alphas (range of composition concentrations)
- Explicit grid/Sobol coverage of the simplex
