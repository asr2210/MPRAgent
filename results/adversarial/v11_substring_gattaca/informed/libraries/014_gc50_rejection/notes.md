# E14: gc50_rejection

## Design
IID random with rejection sampling to per-seq GC ∈ {99, 100, 101}.
~10% acceptance over a 500k pool. Unlike E7 (which flipped bases),
this preserves IID per-position structure exactly.

## Result (mean across 14 evals ≈ 0.7929) — BIG HURT, surprising
- eval_01: 0.8114 (vs E2 0.8565, **-0.045**)
- eval_07 SK-N-SH: **0.3797** (vs E2 0.7723, -0.39) — CRASH
- eval_13: 0.7300 (vs E2 0.8264, -0.10)
- HepG2 ~ 0.90 (preserved)
- K562 ~ 0.83 (preserved)
- SK-N-SH systematically degraded across many evals

## Interpretation — major update to theory
**Tight per-seq GC rejection hurts about as badly as planting motifs.**
This is unexpected: I thought rejection would be a clean test (per-position
marginals preserved). It is not.

**The eval distribution has natural binomial variance in per-sequence GC**.
Random_uniform has GC ~ Binomial(200, 0.5), σ ≈ 7. The eval expects
this spread. Removing the spread (forcing GC ∈ {99,100,101}) is itself
a distribution shift, even though per-position p=0.5 is preserved.

The eval's training distribution is NOT "50% GC sequences"; it's
"sequences drawn from IID Uniform(ACGT)" — which has *natural* GC
variance as a feature.

## Theory v8 — significant refinement
The eval distribution is **IID uniform per position, INDEPENDENT across
positions**. This means:
- Per-position marginals must be uniform.
- Per-sequence GC must follow Binomial(200, 0.5) — NOT a delta.
- Per-sequence base counts must follow the same.
- Any joint constraint (per-seq GC, k-mer diversity, motifs, biology)
  reduces likelihood under IID and hurts.

Random_uniform satisfies all three perfectly. **Random_uniform sits not
just at a local max but at the GLOBAL max of a uniquely-determined
distribution**: IID Uniform^200.

SK-N-SH eval_07 is hypersensitive to ANY joint constraint (it crashes on
gc50_strict, mixed_gc_variance, balanced_bases, rejection, biology).
SK-N-SH is the "joint constraint detector".

## Cross-cell-type implications
Strong implication: this benchmark's eval is mechanically equivalent to
"score a sequence under an oracle that itself was trained on IID Uniform
sequences." Cross-cell-type generalization would mean cross-cell-type
versions of the same oracle. A library optimized for IID would be optimal
across cell types of this same eval family.

For BIOLOGICAL cross-cell-type generalization, this library would lose
badly.

## Next
E15: confirm theory v8 with loose rejection — GC ∈ [90, 110] (about
67% acceptance, captures 1-σ band). Prediction: if v8 is correct, this
should hurt less than E14 (because more natural variance preserved) but
still hurt some. Pinpoints how sensitive the eval is to GC-distribution
tightness.
