# E4: mixed_gc_variance

## Design
- 25,000 IID random uniform (P(N)=0.25 each)
- 25,000 with per-seq GC drawn uniformly from [0.35, 0.65]
- Shuffled, seed 0

## Result (mean across 14 evals ≈ 0.722)
- eval_01: 0.7378 (vs E2 0.8565, **-0.119**)
- All cell types dropped; even K562 and HepG2 lost ~0.05+
- SK-N-SH still dropped on the hardest evals

## Interpretation
Even MODERATE GC variance in HALF the library is severely off-distribution
for the eval. The drop is huge and across all cell types, not just SK-N-SH.

Comparing variance vs performance:
- E3 (per-base balance, zero variance): 0.797
- gc_50 published (zero per-seq GC variance, G:C free): 0.859
- E2 random_uniform (binomial GC variance ~7%): 0.841 (single seed) / 0.857 published
- E4 (mixed wider GC variance ~10-15%): 0.722
- dirichlet (very wide variance): 0.675

Random_uniform sits at or very near the eval distribution. Both narrower
(E3, in per-base sense) and wider (E4) hurt. Published gc_50 marginally
beats random_uniform — confirms per-seq GC fixed at 50% with G:C variance
allowed is the apparent sweet spot.

## Theory v3 → v4
The eval distribution is **statistically very close to IID uniform random
at exactly 50% GC**. Any departure — tighter, wider, structural injection —
causes covariate shift that hurts. The optimum library is approximately
gc_50: per-seq GC fixed at 50% but per-base ratios (G:C, A:T) left free.

E3's failure (G=C strict, SK-N-SH crash) tells us the eval *DOES* sample
G:C imbalance per seq — i.e., the eval distribution is exactly what you'd
get from "draw 100 random GC positions and 100 random AT positions, then
randomly assign each."

**Implication for cross-cell-type generalization**: SK-N-SH's vulnerability
to strict balance suggests it samples sequences with stronger G:C
asymmetry (perhaps neural enhancers tending G-rich or C-rich). The lesson
isn't "use SK-N-SH-specific structure" but "preserve the natural per-position
randomness across G/C and A/T ratios."

## Next
E5 will replicate published gc_50 exactly to (a) verify the small published
lift over random_uniform holds in my pipeline (~+0.003), and (b) set up
a strong control for future small-perturbation tests.
