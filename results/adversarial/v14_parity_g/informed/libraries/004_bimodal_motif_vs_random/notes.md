# 004_bimodal_motif_vs_random

## Design
25,000 motif-packed sequences (4-8 motifs each) interleaved with 25,000
random uniform 200bp sequences. Maximum motif vs no-motif contrast.

## Result
eval_01 mean_r = 0.0026. Still in the noise floor.

## Interpretation
Bimodal label structure doesn't unlock signal. The model cannot reliably
distinguish motif-rich from motif-empty sequences under v14's training
budget, OR motif presence isn't what eval cares about.

Four experiments in (001-004), all ~0. Noise floor confirmed at roughly
[-0.005, +0.005]. Next: pivot to (a) verifying noise floor isn't avoidable
by trying genome-wide random windows and DHS-direct, then (b) committing
to a principled high-diversity library design even if scores stay at floor.
