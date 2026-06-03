# Experiment 008 — Malinois cell-type-discriminating selection

## Design
500k GC-50 random ACGT candidates scored with Malinois. Composite score =
rank_norm(max activity) + 2 * rank_norm(std across 3 cells). Top 50k selected
— prioritizes sequences that activate strongly in some cells but not others.

Selection stats:
- Selected per-seq cross-cell variance: 0.679 (vs random pool ~0.20)
- K562 mean 2.03, std 1.20 (vs 007 magnitude: mean 2.60)
- HepG2 mean 1.75, std 1.05
- SK-N-SH mean 1.90, std 1.42
- GC: mean 0.506, std 0.034

## Hypothesis
Cell-type-DISCRIMINATING sequences expose the model to motif → cell-type
correspondences. This grammar should transfer to held-out cell types (TFs
are shared; only their expression patterns differ).

## Result
- eval_01 = **0.3936** (vs 007 mag 0.3969, 006 span 0.3964) — slightly WORSE
- Mean across 14 evals = **0.3835** (vs 007 = 0.3861, 006 = 0.3860)
- Per-cell-type on eval_01: K562 0.6160, HepG2 0.4303, SKNSH 0.1345
- Runtime: 1235s

## Interpretation
**Cell-type discrimination is NOT the lever.** Selecting for cross-cell
variance actually hurt slightly. The model benefits more from sequences
that are uniformly active across cells (007) than sequences with cell-
type-specific activity patterns.

This is counter to the prior theory that "specificity teaches transfer."
The actual lever is "exposure to high-activity examples," whether they're
broadly active or cell-type-specific.

## Theory v5 update
The Malinois-oracle effect is robustly saturated at ~0.396 regardless of
selection strategy:
| library                | eval_01 | mean_r | strategy |
|------------------------|---------|--------|----------|
| 006 span 5x5x5         | 0.3964  | 0.3860 | uniform 3D coverage |
| 007 top-magnitude      | 0.3969  | 0.3861 | top max(K562,HepG2,SKNSH) |
| 008 cell-type discrim  | 0.3936  | 0.3835 | rank(mag) + 2*rank(std) |

The ~0.005-0.007 lift from oracle is the same regardless of strategy. The
model's representational capacity for "what makes a sequence active" is
quickly saturated at the 50k scale.

## What to try next
- Combine natural sequence statistics with oracle (cCREs + Malinois top)
- Hyperactive generation (simulated annealing / gradient ascent to optimize
  sequences for high Malinois prediction — could give sequences with
  log2FC >> 8 that no random sample can produce)
- Multi-oracle ensemble (Malinois + another prior to disagree)
- Focus on lifting SK-N-SH (still stuck at 0.13-0.15 across ALL libraries)
