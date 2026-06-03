# Experiment 002 — motif-implanted GC-50 random

## Design
50k 200bp sequences: each starts as i.i.d. random ACGT (≈50% GC) scaffold,
then 3 motif instances sampled from random JASPAR PFMs are implanted at
non-overlapping random positions. JASPAR motif median length 9 → ~13% of each
sequence is motif content, rest is random.

## Hypothesis
gc_50 random baseline is 0.397. cCRE balanced (exp 001) was 0.392. If the
model can learn motif content when given a composition-matched scaffold, this
should score above gc_50.

## Result
- eval_01 = **0.3939** (vs gc_50 baseline 0.3972, exp 001 cCRE 0.3919)
- Mean across 14 evals: **0.3827** (vs exp 001 mean 0.3812)
- Per-cell-type on eval_01: K562 0.614, HepG2 0.431, SK-N-SH 0.136
- Runtime: 1217s

## Interpretation
**Essentially identical to experiments 001 (cCRE) and the random baselines.**
Motif content in a clean scaffold did NOT increase learnable signal.

Three libraries with VERY different content all land at 0.39:
| library                     | eval_01 |
|-----------------------------|---------|
| gc_50 random                | 0.397   |
| 001_ccre_balanced (biology) | 0.392   |
| 002_motif_implanted_gc50    | 0.394   |

This strongly suggests the model is NOT primarily learning from motif/biological
features at this training scale. The lever for improvement must be something
else — possibly base composition, activity-range coverage, or specific
sequence-feature distributions that match the eval set.

## Per-cell-type pattern (identical across libraries!)
- K562: 0.60 - 0.61
- HepG2: 0.43 - 0.43
- SK-N-SH: 0.14 - 0.15

The cell-type bias is a TRAINING-INVARIANT phenomenon for these strategies.
The model converges to similar per-cell-type performance regardless of what
training sequences are provided. Likely the eval set difficulty per cell type
is the dominant factor.

## Next step
Need to test something more dramatic. Possible directions:
- Pure promoter (PLS-only) library — promoters drive HIGH activity in MPRA
- Activity-range stratified library — explicitly span quiet→active
- HepG2/SKNSH-specific elements — lift the lagging cell types
- Mix of strategies with VERY different distributions
