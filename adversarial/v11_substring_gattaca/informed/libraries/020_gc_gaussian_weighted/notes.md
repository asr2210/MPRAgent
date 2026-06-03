# E20: gc_gaussian_weighted

## Design
IID random with importance-weighted acceptance: accept each with
probability exp(-(GC-100)²/(2·5²)). Soft Gaussian taper around GC=100,
σ=5. Vs E15's hard cutoff [90, 110].

## Result
- eval_01: 0.8640 (vs E15 0.8777, -0.014)
- mean_r: 0.8436 (vs E15 0.8591, -0.015)
- SK-N-SH eval_07: 0.5914 (vs E15 0.7183, -0.13)
- Worse across most evals

## Interpretation
Gaussian weighting is more peaked than uniform-on-[90,110]. The effective
SD is ~5.2 (similar to E15 5.23) but the weight distribution is peaked at
100. SK-N-SH crashes again, indicating too much concentration.

**Hard cutoff at [90, 110] beats Gaussian taper**. The eval likely has
hard or near-uniform GC band, not Gaussian. (Or the model's preferred
training distribution is uniform on the band — unclear why exactly, but
hard cut wins.)

## Next
E21: layer RC pairs on GC band. 25k GC-filtered + their 25k RCs.
Tests if RC augmentation (which preserves GC) adds value when combined
with GC band.
