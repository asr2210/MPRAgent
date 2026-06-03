# E21: gc_band_rc_pairs

## Design
25k IID random with GC ∈ [90, 110] + their 25k reverse complements
(RC preserves GC count). Layered RC augmentation on E15's optimum.

## Result
- eval_01: 0.8759 (vs E15 0.8777, -0.002)
- mean_r: 0.8570 (vs E15 0.8591, -0.002)
- Essentially matches E15 within noise

## Interpretation
RC pairing adds nothing on top of GC band. E8 already showed RC matches
random; now confirmed it also matches GC-filtered random. The model
likely has built-in RC equivariance, so explicit RC pairs are redundant.

## Next
E22: variance check — E15's design with seed=1. Confirms reproducibility
of the +0.02 gain (vs single-seed noise ~0.006).
