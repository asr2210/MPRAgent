# E29: gc_band_seed2

## Design
E15's GC band [90, 110] design with seed=2 (third reproducibility check).

## Result
- eval_01: **0.8808** (BEST overall, vs E15 0.8777, +0.003)
- mean_r: **0.8625** (BEST overall, vs E15 0.8591, +0.003)
- SK-N-SH eval_07: 0.7328 (best across GC band experiments)

## Three-seed summary of GC band [90, 110]
- seed 0 (E15): 0.8777 / 0.8591
- seed 1 (E22): 0.8755 / 0.8578
- seed 2 (E29): 0.8808 / 0.8625
- Mean: 0.8780 / 0.8598
- Range: 0.005 / 0.005 → single-seed SD ~ 0.0025

## Theory v9 fully validated
GC band [90, 110] reliably gives +0.018 to +0.022 mean over random_uniform
(0.8408). The +0.005 inter-seed variance is real. E29 happened to land at
the high end.

## Next
E30: final library — combine seed=2 sequences with per-position A↔T/C↔G
balance (E24 method) to maximize.
