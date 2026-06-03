# 017 — Full JASPAR + cell-type-weighted draws (60/40)

## Plan
Load all 2344 JASPAR vertebrate PFMs. Per insertion, 60% draw from
cell-type-targeted subset (289 PFMs), 40% from full pool. Same 3 motifs/seq
random uniform backbone otherwise.

## Result
**eval_01 = 0.4195.** Within noise of 008 plateau (-0.009).
K562: 0.585 (↓ 0.011), HepG2: 0.616 (↓ 0.013), **SK-N-SH: 0.057** (↑ 0.007).
eval_07 SK-N-SH = **0.069**, eval_10 = 0.060, eval_13 = 0.065 (multiple
evals breaking 0.06).

## What this teaches
- **Trade-off pattern:** broader pool hurts K562/HepG2 specificity but
  benefits SK-N-SH. Net effect on mean_r is roughly flat.
- Confirms T11: SK-N-SH oracle responds to broader sequence variety, not
  just neural TFs.
- "SK-N-SH floor at 0.06" is NOT a structural property of the oracle. It
  can be moved by library composition. Highest single-eval seen so far:
  0.069 (eval_07, this experiment).

## Theory T13
The three cell types have DIFFERENT preferences:
- K562/HepG2 want narrower, cell-type-targeted motif pools (sharper signal).
- SK-N-SH wants broader composition (variety across the sequence pool).

These preferences pull in opposite directions, capping mean_r near 0.42-0.43
under a single-strategy library. To break the plateau, may need a STRUCTURED
library that gives each cell type what it wants — e.g., some sequences with
K562-only sharp motifs, others with neural-broad inserts.

## Generalization implication
The "broader pool helps generalization" intuition is half right:
broader → better for noisy/composition-sensitive oracles (SK-N-SH-like
unknowns), worse for sharp/motif-sensitive ones (K562/HepG2-like
unknowns). For unknown cell types, library should hedge: include both
narrow and broad sequences.

## Next
Exp 018: STRUCTURED per-sequence balance — every sequence contains EXACTLY
1 K562 motif + 1 HepG2 motif + 1 SK-N-SH motif (one from each cell-type
pool). Tests whether forced per-seq cell-type balance beats random uniform
draws from the combined pool. Hypothesis: gives each cell-type oracle a
clear signal in every training sequence.
