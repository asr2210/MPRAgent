# Experiment 010 — Gosai UKBB-only

## Method
50K from UKBB-only Gosai (data_project=="UKBB"), lfcSE<0.3, quintile stratified.
255K UKBB sequences passed quality filter.

## Result
**eval_01 = 0.0130.** Below the 0.018 plateau.

## Interpretation
UKBB-only underperforms mixed UKBB+GTEX. Either:
- GTEX contributes useful signal (different variant types / genomic contexts)
- Library size of 255K is still enough to sample from, so it's composition
- UKBB is GWAS-disease variants only, narrower distribution than GTEX

K562 is the strongest cell type here (0.022), suggesting UKBB is K562-heavy.

Next: GTEX-only as a parallel ablation.
