# E19: gc_and_perbase_band

## Design
IID random rejection to per-seq GC ∈ [90, 110] AND each base count ∈
[40, 60]. Adds per-base filter on top of GC.

## Result
- eval_01: 0.8757 (vs E15 0.8777, -0.002)
- mean_r: 0.8564 (vs E15 0.8591, -0.003)
- SK-N-SH eval_07: 0.6744 (vs E15 0.7183, -0.04)
- Slight loss across most evals

## Interpretation — per-base filter doesn't help on top of GC
Per-base filtering is a stricter constraint that removes sequences the
eval distribution might still like. GC alone is the right macro feature
to filter on.

This is informative: the eval doesn't filter on per-base count
independently of GC. It's a clean GC-only soft filter.

## Next
E20: try Gaussian-weighted GC (soft taper) instead of hard cutoff.
If the eval distribution has a SOFT GC band, Gaussian weighting will
match it better than uniform-on-[90,110]. If hard cutoff is already
optimal, Gaussian won't help.
