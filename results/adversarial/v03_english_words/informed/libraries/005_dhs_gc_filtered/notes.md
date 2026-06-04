# 005 — DHS filtered to 45-55% GC

## Plan
Test whether DHS underperformance vs random uniform is due to GC variance.
Filter to GC ∈ [0.45, 0.55] (matches random_uniform's natural range).

## Result
**eval_01 = 0.3876.** Indistinguishable from unfiltered DHS (0.3868).

Library GC: 0.493 ± 0.030 (tight, as designed). Still 0.39 eval_01.

## Conclusion
GC variance is NOT the cause of DHS underperformance. Whatever the pipeline
penalizes about DHS sequences is something else — perhaps:
- Repeat content (Alu, LINE, low-complexity regions natural in DHSs)
- Specific dinucleotide composition (CpG islands, etc.)
- k-mer distribution mismatch with the surrogate's training data

OR more fundamentally: the pipeline's surrogate model is too simple to extract
value from natural sequence structure, and treats DHS sequences as "noisy
random" inputs that hurt training stability.

## Theory T5
Strategies.md baselines suggest a CEILING around 0.43 for "simple structured"
libraries (gc_50 best at 0.4243). DHS libraries hit ~0.39 (worse). This
implies the pipeline rewards a specific distribution that random uniform
approximates but DHS does not. The next informative test is whether any
*content-rich* synthetic library (random backbone + injected regulatory
motifs) can break the 0.43 ceiling.

## Next
Exp 006: Random uniform backbone with JASPAR motif insertions. If motif
content adds informative variance, this should beat random_uniform (0.42)
and break the 0.43 strategies.md ceiling. If not, the pipeline is at hard
ceiling for biological content.
