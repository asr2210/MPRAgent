# E16: gc_band_wider

## Design
IID random rejection to per-seq GC ∈ [85, 115] (~86% acceptance band).
Wider than E15's [90, 110]. Tests if optimum is wider or narrower than 1-σ.

## Result
- eval_01: **0.8688** (vs E15 0.8777, -0.009; vs E2 0.8565, **+0.012**)
- mean_r: 0.8451 (vs E15 0.8505, -0.005; vs E2 0.8408, +0.004)
- All evals improved over E2 but slightly worse than E15

## Interpretation — sweep is finding an optimum
GC band sweep so far:
- [99, 101]  → 0.7929 mean  (E14, too tight)
- [90, 110]  → **0.8505 mean** (E15, peak so far)
- [85, 115]  → 0.8451 mean  (E16, looser, slightly worse)
- (no band)  → 0.8408 mean  (E2 random, baseline)

The optimum is at or tighter than [90, 110].

## Next
E17: GC ∈ [93, 107] (tighter than E15). If improves further, peak is
narrower than 1-σ. If degrades, peak is in [90, 110] neighborhood.
