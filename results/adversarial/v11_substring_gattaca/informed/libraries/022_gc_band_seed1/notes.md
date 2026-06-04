# E22: gc_band_seed1

## Design
E15's design (GC ∈ [90, 110] rejection) with seed=1.

## Result
- eval_01: 0.8755 (vs E15 0.8777, -0.002, within noise)
- mean_r: 0.8578 (vs E15 0.8591, -0.001)

## Interpretation
E15's gain is **reproducible across seeds**. Single-seed noise is ~0.006,
so this match within 0.002 confirms the GC-band optimum.

## SK-N-SH eval_07 observation
SKNSH eval_07 in GC-band experiments runs 0.69-0.74 — slightly LOWER than
random_uniform (0.77). The GC band hurts SK-N-SH eval_07 specifically
while helping all other evals enough to net +0.018 mean.

Hypothesis: SK-N-SH eval_07 prefers sequences with broader GC variance
(including extremes). Could test by mixing GC-band + full random.

## Next
E23: 25k GC-band + 25k full random mix. Tests if SK-N-SH eval_07
recovers while keeping most of the GC-band improvement.
