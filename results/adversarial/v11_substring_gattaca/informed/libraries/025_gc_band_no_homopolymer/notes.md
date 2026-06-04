# E25: gc_band_no_homopolymer

## Design
E15's GC band [90, 110] + reject any sequence with 10-bp homopolymer.
Filter rejected only ~50 of 86k post-GC sequences (~0.06%).

## Result
- eval_01: 0.8775 (vs E15 0.8777, essentially tied)
- mean_r: ~0.859 (same as E15)
- No effect; homopolymer outliers are too rare to matter

## Interpretation
Local-bias outliers don't affect training. E15 is robust.

## Next
E26: per-position EXACT balance only (no GC filter). Tests whether the
gain in E24 came from GC filter or per-position balance.
