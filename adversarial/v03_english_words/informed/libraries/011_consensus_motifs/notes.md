# 011 — Cell-type-targeted CONSENSUS motifs (no PFM stochasticity)

## Plan
Same as exp 008 (best, 0.4283): 289 cell-type-targeted PFMs, 3 motifs per 200bp
random uniform backbone. Difference: insert deterministic consensus string
(argmax at each position) instead of stochastically sampled PFM instance.
Tests whether PFM noise was diluting signal.

## Result
**eval_01 = 0.4169.** Worse than exp 008 (0.4283) by -0.011.
K562: 0.588 (↓ from 0.596), HepG2: 0.618 (↓ from 0.629), SK-N-SH: 0.045 (≈).

## What this teaches
- **PFM stochasticity HELPS, not hurts.** Pure consensus repeats the same
  string for the same TF every time. Stochastic sampling produces diverse
  motif instances (different positions of the consensus, occasional flanking
  variations), which gives the surrogate more training samples per TF.
- Consensus inserts effectively reduce vocabulary size — 289 unique 8-20mer
  strings vs ~tens of thousands of distinct stochastic samples. The
  surrogate over-fits the exact consensus and generalizes less to held-out
  motif instances.
- Motif identity matters, but motif INSTANCE DIVERSITY matters too.

## Theory T9
T8 said information density per sequence matters. T9 adds: per-motif
INSTANCE diversity also matters. A library should expose the surrogate to:
1. Many cell-type-relevant TFs (signal)
2. Many distinct instances of each (generalization)
PFM-sampling does both naturally; consensus-only collapses (2).

## Next
Exp 012: same 289-motif pool but force multi-INSTANCE of the same TF in a
single sequence (e.g., 3x stochastic samples of GATA1 in one seq). Tests
whether clustering of same-TF instances within a sequence creates stronger
per-TF features. If this helps: the surrogate learns motif "density" per
TF, not just presence/absence.
