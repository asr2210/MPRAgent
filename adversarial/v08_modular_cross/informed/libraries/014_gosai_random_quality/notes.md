# Experiment 014 — Random Gosai with quality filter (no stratification)

## Method
Random 50K from Gosai lfcSE<0.5 pool (697K sequences). No stratification.

## Result
**eval_01 = 0.0168** (vs 0.0181 stratified, vs 0.0000 random no-filter).

## Interpretation — KEY FINDING
Stratification adds essentially nothing on top of quality filter.
- Quality filter explains the lift from 0.000 → 0.017
- Stratification adds ~0.001 (within noise)

So the dominant axis is QUALITY (lfcSE filter), not activity distribution shape.
This explains why every selection-on-activity variant (extremes, variance,
percentile shifts) underperforms — they're trading off against the quality
filter's natural sample.

Also: this library has the highest eval_07/08/13 values for Gosai libraries
so far (0.005/0.005/0.010) — more diversity helps the "hard" eval cluster.

## Theory
The model needs CLEAN LABELS on a NATURALLY DISTRIBUTED training set. Quality
filter removes noisy labels; stratification distorts the distribution. The eval
distribution must roughly match the post-quality-filter Gosai distribution.

To break the ceiling: improve quality further without losing diversity. Options:
- max-SE filter (stricter on all 3 cells)
- Multiple lfcSE thresholds per cell (per-cell filtering)
- Data augmentation (revcomp, shuffles)
- External substrate combined with Gosai
