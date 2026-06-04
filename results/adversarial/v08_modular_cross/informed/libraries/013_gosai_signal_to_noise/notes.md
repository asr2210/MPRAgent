# Experiment 013 — Signal-to-noise selection

## Method
Score each Gosai sequence by mean |log2FC| / (lfcSE + 0.05) across 3 cell types
(with lfcSE<0.5 filter). Take top 100K by SNR, quintile-stratify by mean
activity to keep diversity. Cuts noise without losing signal.

## Result
**eval_01 = 0.0126** — below plateau.

## Interpretation
High-SNR sequences are essentially "loud and clean". They overlap heavily with
activity extremes (since |effect| is in the numerator). This recapitulates the
extremes-selection bias: boosts eval_04/09 (0.0185) but loses eval_01.

## Pattern across experiments
Every selection axis I've tested (activity extremes, cell-type variance,
sub-source isolation, tight quality, SNR) underperforms broad stratification on
eval_01. The 0.018 plateau is firm for single-criterion Gosai libraries.

## Next direction
Two unexplored axes:
1. Activity-distribution MATCHING (not just stratifying) — match a natural,
   non-uniform shape that may better mirror the eval's own composition
2. Off-Gosai augmentation — REF/ALT pairs, k-mer diversity, secondary MPRA
