# 006 — Random + 3 JASPAR motifs/sequence

## Plan
Test: does adding JASPAR TF motifs to random uniform backbones beat plain
random uniform (0.42)? Each sequence: random uniform 200bp + 3 motifs sampled
from JASPAR2024 (2344 vertebrate PFMs) inserted at random positions.

## Result
**eval_01 = 0.4252.** Marginally above random_uniform (0.4203) by ~1.2% absolute.
Consistent across all 14 evals. K562 went 0.585 → 0.591, HepG2 0.617 → 0.622,
SK-N-SH stuck at 0.06.

## What this means
- The surrogate CAN exploit motif content — motifs are "learnable" features.
- Motifs marginally help K562/HepG2 but NOT SK-N-SH (which is stuck regardless).
- 0.42 was NOT a hard pipeline ceiling — biology can break it.

## Implications for theory
T6: Motif content provides informative variance for K562/HepG2 prediction
(small effect size, ~1% absolute). SK-N-SH is a hard floor problem.

The improvement is small because:
- 3 motifs in 200bp is a low density (~5% of sequence is motif)
- Random selection from 2344 motifs dilutes signal
- The surrogate might learn shallow features and saturate quickly

## Next
Two paths:
A) **Push motif density**: 5-8 motifs per sequence
B) **Cell-type-specific motifs**: focus on K562/HepG2/SK-N-SH-relevant TFs

Let me try both — exp 007 tests density (6 motifs/seq), exp 008 tests
cell-type-targeted motifs.

Also: switch backbone from random uniform to gc_50 (50% GC, base-shuffled).
gc_50 baseline (0.4243) was slightly above random_uniform (0.4228), so gc_50
backbone + motifs should be a slightly better foundation.
