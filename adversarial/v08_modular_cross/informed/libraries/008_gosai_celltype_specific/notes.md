# Experiment 008 — cell-type-specific (high variance) Gosai

## Method
50K Gosai with lfcSE<0.3, ranked by variance across (K562, HepG2, SKNSH).
Take top 50K by cell-type variance.

## Result
**eval_01 = 0.0122** (vs 005 plateau 0.0180). Lower than stratified.

## Interpretation
Cell-type-specific selection doesn't help. The eval distribution is not
primarily driven by cell-type-discriminating sequences. Broad activity
coverage (stratification) wins.

## Theory update
Specialized selections (extremes, cell-type-specific) all underperform broad
stratification for eval_01. The plateau at ~0.018 appears to be the ceiling
for "quality-filtered Gosai with reasonable coverage." Pushing past requires
either:
- Different data source (other MPRA datasets, synthetic)
- Larger N (locked at 50K)
- More informative per-example signal
