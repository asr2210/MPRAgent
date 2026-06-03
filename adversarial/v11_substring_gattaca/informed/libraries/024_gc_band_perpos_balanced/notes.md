# E24: gc_band_perpos_balanced

## Design
E15's GC band [90, 110] pool, then per-position rebalance via A↔T and
C↔G swaps (preserves per-seq GC). Reduces per-position sampling noise.

## Result
- eval_01: 0.8774 (vs E15 0.8777)
- mean_r: 0.8587 (vs E15 0.8591)
- Effectively tied with E15

## Interpretation
Per-position sampling noise (~0.78% deviation) is too small to affect
performance. The natural variance is already at the noise floor.

## Theory v9 holds robustly
GC band [90, 110] is the ONE filter that matters. All other tested
filters (per-base, k-mer entropy, RC pairs, Gaussian, balanced) add
nothing.

## Next
E25: GC band + reject sequences with ≥10-bp homopolymer runs. Tests if
removing local-bias outliers helps further or hurts.
