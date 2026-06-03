# Experiment 003 — PLS-only library

## Design
50k 200bp windows from ENCODE cCRE PLS (Promoter-Like Signature) elements only.
~41k PLS in cCRE pool; extras generated via random jitter ±50bp from element
center.

## Hypothesis
Promoters drive the strongest MPRA signal. A promoter-focused library should
provide a richer activity training signal than balanced/random sequences.

## Result
- eval_01 = **0.3748** — WORSE than cCRE balanced (0.3919) and gc_50 (0.397)
- Mean across 14 evals: **0.3653**
- Per-cell-type: K562 0.58, HepG2 0.41, SK-N-SH 0.13
- Runtime: 1210s

## Why it's worse: GC composition penalty
The PLS library has **mean GC = 0.617**, far from the optimum 0.50. CpG island
density in promoters drives this. Strategies.md shows:
- gc_50 (50% GC): 0.397
- gc_rich (80% GC): 0.222
PLS at 62% GC sits between, scoring 0.375 — consistent with composition cost.

## Big update to theory
**Composition is the DOMINANT lever** at this 50k-library scale on this evaluator.
The model performance closely tracks training-set base composition:
| library                  | mean GC | eval_01 |
|--------------------------|---------|---------|
| gc_50 (baseline)         | 0.50    | 0.397   |
| 002 motif-implanted gc50 | 0.50    | 0.394   |
| 001 cCRE balanced        | 0.51    | 0.392   |
| 003 PLS only             | 0.62    | 0.375   |
| at_rich (baseline)       | 0.20    | 0.299   |
| gc_rich (baseline)       | 0.80    | 0.222   |

Performance drops as composition deviates from 50%. Biological content provides
no measurable lift above what composition predicts.

## Next experiment idea
Test whether GC-FILTERED cCREs (selecting only those with GC ∈ [0.45, 0.55])
can BEAT gc_50. If they can, biology adds value when composition is matched.
If not, composition is the only lever and biology truly doesn't help here.
