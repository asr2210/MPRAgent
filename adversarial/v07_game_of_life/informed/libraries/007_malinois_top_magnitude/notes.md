# Experiment 007 — Malinois TOP-magnitude (high-activity only)

## Design
500k GC-50 random ACGT candidates scored with Malinois. Selected top 50k by
max(predicted log2FC) across K562, HepG2, SK-N-SH. No span control —
exclusively the highest-activity tail.

Selection stats:
- Selected mean predicted activity: K562 2.60, HepG2 2.36, SKNSH 2.63
  (vs full pool means: K562 0.70, HepG2 0.65, SKNSH 0.67)
- Min selected max-cell score: 2.32; max: 9.66
- Selected GC: mean 0.511, std 0.034

## Hypothesis
Tests if exp 006's lift came from SPAN (diverse activities) or MAGNITUDE
(presence of high-activity sequences). If 007 > 006, magnitude dominates.

## Result
- eval_01 = **0.3969** (vs 006 = 0.3964, gc_50 = 0.397)
- Mean across 14 evals = **0.3861** (vs 006 = 0.3860)
- Per-cell-type: K562 0.618, HepG2 0.436, SKNSH 0.137
- Runtime: 1228s

## Interpretation
**Tied with exp 006.** Span and top-magnitude give identical mean_r and very
similar per-cell-type scores. This means:
- The +0.005 lift from oracle selection is robust
- It saturates quickly — both extremes (span+magnitude) give the same answer
- The lift comes from EXPOSURE TO HIGH-ACTIVITY EXAMPLES, not the diversity
- Adding low-activity sequences does NOT help (top-only ≈ span which includes
  low-activity)

This implies the bottleneck is the model's ability to predict the
high-activity tail. Once trained on enough high-activity sequences, the
predictions stabilize.

## Theory v5 update
The Malinois-oracle-selection lever is real but saturated at ~0.396-0.397.
Both span and magnitude reach the same plateau, slightly above the
composition-only baseline. To break further:

1. Cell-type DISCRIMINATION — select for high variance across cells
2. Natural sequence + oracle — see if cCRE-derived top-activity > random-derived
3. Hyperactive design via FastSeqProp / SimulatedAnnealing — generate
   sequences explicitly optimized for high predicted activity
4. SK-N-SH-specific motif augmentation — only lever yet unexplored for SKNSH
