# E28: gc_uniform_on_band

## Design
Per-seq GC distribution UNIFORM on [90, 110] (~2381 sequences at each GC
value 90..110). vs E15's truncated-Binomial shape (peaked at GC=100).

## Result
- eval_01: 0.8713 (vs E15 0.8777, -0.006)
- mean_r: 0.8540 (vs E15 0.8591, -0.005)

## Interpretation
Uniform-on-band loses slightly to truncated-Binomial. The eval prefers
the naturally-peaked distribution (more mass near GC=100). E15's
rejection sampling — which preserves the natural shape within the band —
is optimal for the eval.

## Next
E29: third variance check (seed=2) on E15's GC band [90, 110] design.
Confirms reproducibility (E15 0.8777, E22 0.8755, expect E29 ~0.876).
