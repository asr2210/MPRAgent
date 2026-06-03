# E18: gc_band_refine

## Design
IID random rejection to per-seq GC ∈ [88, 112]. Bisects between
E15 [90,110] and E16 [85,115].

## Result
- eval_01: 0.8757 (vs E15 0.8777, -0.002, virtually tied)
- mean_r: 0.8579 (vs E15 0.8591, -0.001)
- SK-N-SH eval_07: 0.7379 (vs E15 0.7183, +0.02)

## Interpretation — flat plateau confirmed
GC band optimum is a flat plateau between [88, 112] and [90, 110].
Both give ~0.858 mean / ~0.876 eval_01. Differences within noise.

Sweep wrap:
- [99, 101]  → 0.79  (E14)
- [93, 107]  → 0.85  (E17, descending)
- [90, 110]  → 0.859 (E15, peak)
- [88, 112]  → 0.858 (E18, plateau)
- [85, 115]  → 0.852 (E16, mild descent on wide side)
- random     → 0.841 (E2)

**GC band optimum locked in at [90, 110]**. Move to layered improvements.

## Next
E19: layer per-base count constraint on top of GC band. Filter to
GC ∈ [90, 110] AND each base count ∈ [40, 60]. Tests whether the eval
has additional per-base filtering beyond GC, or whether GC subsumes it.
