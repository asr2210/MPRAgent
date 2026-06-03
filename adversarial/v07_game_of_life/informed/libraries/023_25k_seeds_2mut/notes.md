# Experiment 023 — 25k unique random + 1 mutant each (independence bridge)

## Design
25k unique GC-50 random + 25k 1-base mutants of those = 50k.
Each pair (seed, mutant) shares 199/200 positions.

## Result
- eval_01 = **0.3902** (Δ vs 007 = -0.0067)
- mean_r = **0.3799** (Δ vs 007 = -0.0062)
- Per-cell: K562=0.615, HepG2=0.427, SK-N-SH=0.128

## Interpretation
Slight drop vs 007 (~0.006), borderline of noise. Compare:
| Config | Effective indep | mean_r |
|--------|----------------|--------|
| 50k indep random (007) | 50k | 0.3861 |
| 25k indep + 25k RC (020) | 25k (RC-aug) | 0.3861 |
| 25k indep + 25k 1-base mut (023) | ~25-30k | 0.3799 |
| 5k indep + 9 1-base mut (021) | ~5k | 0.3600 |

**Resolution:** RC pairs are handled by trainer's built-in RC augmentation
(so effectively 1 sample per pair, but model learns RC-equivariance "for free").
1-base mutants are treated as NEW data — the trainer expends gradient steps
on them but they provide redundant info, slightly hurting compared to
truly independent samples.

The drop from 50k→25k effectively-independent in 023 (-0.006) is much
smaller than the 50k→5k drop in 021 (-0.026), confirming roughly
logarithmic scaling of info content.

## Theory v15 reinforced
The relationship between effective independent samples and mean_r:
- 5k → 0.360
- 25k → 0.380-0.386 (with various augmentations)
- 50k → 0.386

Diminishing returns saturate around 25k indep samples.
