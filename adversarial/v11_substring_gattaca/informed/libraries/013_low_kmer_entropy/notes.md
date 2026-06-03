# E13: low_kmer_entropy

## Design
Same 200k random pool as E12. Selected BOTTOM 50k by 6-mer diversity
(unique k-mer count 170-189, mean ~187). Symmetric opposite of E12's
top-quartile selection.

## Result (mean across 14 evals ≈ 0.8296)
- eval_01: 0.8449 (vs E2 0.8565, **-0.012**)
- mean_r: 0.8296 (vs E2 0.8408, -0.011; vs E12 0.8428, -0.013)
- HepG2 dropped most (~0.86 vs E12's 0.90)
- SK-N-SH eval_07 ~ 0.69 (same as E12)

## Interpretation — clean asymmetry result
E12 (top diversity) was flat vs random. E13 (bottom diversity) hurts by
~0.011 mean. This is **small but above the ~0.006 noise floor**.

**Implications**:
1. The eval IS sensitive to micro-distributional structure, not just
   per-position marginals.
2. Asymmetry: deviation from random distribution always hurts, never
   helps. Random sits at the local maximum.
3. Selecting low-diversity sequences likely introduces subtle base-bias
   (k-mer repetition correlates with run-length biases), and the model
   picks this up.

## Theory v7 — refined
Eval = synthetic, very near IID random uniform 50% GC. Random_uniform is
at the local maximum of the manifold. **Every direction off random hurts.**
Magnitude of hurt scales with magnitude of deviation:
- Tiny (low-kmer): -0.01
- Modest (planted motifs, mutation-pairs): -0.02
- Larger (Markov, gc50_strict): -0.04 to -0.06
- Major (mixed_gc_variance, real biology): -0.07 to -0.13

The benchmark ceiling for the eval_01 metric is ~0.857. **No tested
deviation has improved on random_uniform — and now we know symmetric
shrinkage from random hurts more in one direction than the other.**

## Cross-cell-type implications
The asymmetry shows the eval isn't measuring "generalization potential"
in a meaningful biological sense — it's measuring "match to random
uniform distribution". A library that better matches biology would lose
on this eval, definitively.

## Next
E14: try rejection-sampled exact GC=100 random. Tests whether the
gc_50 published baseline's tiny +0.0025 edge over random_uniform is
reproducible or noise. Quick, decisive.
