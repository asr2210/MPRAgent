# Skill: Decomposing this MPRA pipeline's eval ceiling

## TL;DR
This `prepare.py` caps at eval_01 ≈ 0.34 regardless of library design.
Per-cell caps: K562 ≤ 0.17, HepG2 ≤ 0.22, SKNSH ≤ 0.63.
The model behaves as a dinucleotide-frequency regressor.

## How to decompose any pipeline's input dependency
Use these four probe libraries to map where the lift comes from:

| Probe                       | Tests                                |
|-----------------------------|--------------------------------------|
| Uniform random {A,C,G,T}    | Baseline / noise floor               |
| GC-distribution-matched random | Lift from mononucleotide composition |
| Mononucleotide-shuffled real | Same (control — should match above) |
| Dinucleotide-shuffled real  | Lift from 2-mer structure            |
| Intact real sequences       | Lift from grammar/motifs             |

If shuffled ≈ intact, the model isn't using motif info. If GC-matched ≈
shuffled, the model only uses GC%.

## Specific numbers found in informed_claude run
- Uniform random:        0.244
- GC-matched random:     0.306 (+0.062)  ← biggest lever
- Dinuc-shuffled Gosai:  0.320 (+0.014)
- Intact Gosai random:   0.323 (+0.003)  ← essentially no grammar lift
- + chr 7/9/13/21/X:     0.336 (+0.013)  ← Gosai-quirk
- + top-SNR within chr:  0.342 (+0.006)

## eval-set twin pairs (informed_claude pipeline)
- eval_01 = eval_05 (the primary metric has a twin)
- eval_02 = eval_05 = eval_14
- eval_03 = eval_12
- eval_04 = eval_09
- eval_06 = eval_11
- eval_08 is broken (~0.08) across all libraries

## Practical implication
Library design effort beyond "match Gosai's GC distribution and use
chr 7/9/13/21/X high-SNR sequences" gives diminishing returns.
Don't spend time on:
- Motif engineering
- Source mixing  
- RC augmentation
- Stratification
- Sequence duplication
