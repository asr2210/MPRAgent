# Experiment 010 — Dirichlet(0.3)

## Result
eval_01 = **0.1344** (vs dirichlet(0.5) 0.1395 — worse by 0.005)

## Interpretation
Alpha=0.3 is too extreme. More near-homopolymer sequences in the tail
hurt performance (consistent with v06 homopolymer_rich = 0.058).

Alpha sweet spot is around 0.5:
| alpha | eval_01 |
|-------|---------|
| 0.3   | 0.1344  |
| 0.5   | 0.1395  |
| 1.0   | 0.1371 (Sobol-uniform = effective alpha ≈ 1) |
| 0.3-2.0 | 0.1377 (mixed) |

The Dirichlet(0.5) optimum is robust.

## Next
Try explicit composition stratification: generate many dirichlet
candidates, ensure equal representation across (dominant, second)
composition modes. Tests if dirichlet(0.5) undersamples some modes.
