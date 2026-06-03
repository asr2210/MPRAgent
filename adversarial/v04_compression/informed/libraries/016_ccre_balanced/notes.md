# 016_ccre_balanced

## Design
12.5k each from 4 ENCODE cCRE element classes:
- PLS (promoter-like)
- pELS (proximal enhancer-like)
- dELS (distal enhancer-like)
- other (CTCF-only / DNase-H3K4me3)
200bp window centered on each element midpoint. seed=0.

## Result — strong negative
- eval_01 = **0.4022** (vs E2 0.4992; Δ -0.097)
- mean over 14 ≈ 0.394
- Lower than E3 NMF-stratified DHS (0.438) and DHS uniform (0.470).

This is the worst real-DNA result so far. Pretty close to motif-implant
(E7 = 0.347) and Markov k=4 (E5 = 0.268) territory.

## Interpretation
Forced equal-class balance OVER-REPRESENTS the rare/narrow classes
(promoter ~4% naturally → 25% here). E9 already showed promoter-only
= 0.37. With 25% promoter weight, the library inherits that penalty.

Element-class stratification is the WORST kind of stratification I've
tested. NMF-stratification (E3) gave 0.438 (4-5x finer-grained
strata but each class is at least DHS-related). cCRE classes are
coarser but include the pathologically narrow promoter class at 25%.

This is consistent with the theory: any stratification that over-
represents a narrow distribution within real DNA hurts.

## Theory update v6.2
The "over-weighting rare narrow classes" penalty is the dominant
cost of stratification. Naturally proportional sampling avoids it.

Summary of stratification penalties:
- Uniform DHS sampling:       0.470 (mild narrowing only)
- NMF-stratified DHS:         0.438 (NMF over-weights rare programs)
- cCRE-class balanced:        0.402 (over-weights rare promoters)
- Promoter only (extreme):    0.370 (the over-weighted class alone)

There is a hierarchical penalty: the more aggressively a curation
over-represents a narrow class, the worse it gets.

## What's left

The 0.50 ceiling is the natural cost of "50k slots from broad human".
To break it, I need to either:
(a) Improve the QUALITY of random slots (filter out wasted slots
    like pure-repeat / N-rich tiles).
(b) Find a different overall distribution closer to the eval set.
(c) Use multi-seed and hope for sample variance (but E14 showed this
    is <0.001).

## Plan for E17
**Quality-filtered random**: 50k tiles from chr1/17/19/22 but reject:
- any tile with a single nucleotide >40% of 200bp (simple repeat / N)
- bottom 10% of dimer entropy (tandem repeats)

Hypothesis: ~20-30% of random tiles are "wasted slots" in low-info
regions. Replacing them lifts score by 0.01-0.03.

Outcomes:
- ≥0.51: quality filter is the path forward; iterate stronger filters.
- ~0.50: neutral; median random is already informative.
- <0.49: filter is too aggressive, hurts diversity.
