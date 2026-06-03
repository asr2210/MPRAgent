# 014_dhs_motif_augmented — notes

## Design
dhs_specific sequences with 4 random JASPAR motifs planted at random
positions (overwriting real bases).

## Result
eval_01 = 0.0430  (vs dhs_specific 0.049) — DROP
mean_r = 0.0402

## Interpretation
Motif augmentation HURTS. Inserting synthetic JASPAR motif samples
disrupts the real motif content already present in DHS sequences.

This combined with Exp 013 (shuffle drops 0.049 → 0.038) means:
The eval rewards REAL biological sequence structure, NOT motif density.
Real DHS sequences already have a "right amount" of motif content; both
adding more synthetic and destroying existing both hurt.

Implication: stay close to real biological sequences. Optimization is
in WHICH biological sequences to pick, not in modifying them.
