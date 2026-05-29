# 024 — 3-source synthetic split

## Goal
Test if splitting 15k synthetic across three sources (motif, uniform random,
homotypic cluster — 5k each) gives compound diversity benefit.

## Method
- 67,500 cCRE-fwd + 67,500 different cCRE-RC (135k unique, mixed strand)
- 5,000 motif-embedded + 5,000 uniform random + 5,000 homotypic clusters
- Total 150k, seed=23

## Result: 0.8861 — within noise of cluster but on the low end.

| design | mean |
|--------|------|
| 015 (single motif synth) | 0.8905 |
| 016 (FANTOM half + motif half synth) | 0.8907 |
| 022 (replicate of 015) | 0.8875 |
| 023 (kitchen sink) | 0.8886 |
| 024 (3-source synth) | 0.8861 |

## Key observations
1. **Synthetic source diversity does NOT compound.** Three-way mix lands
   at low end of the noise cluster.
2. **Possible weak effect**: splitting synthetic into smaller chunks
   may slightly dilute each chunk's contribution (each at 5k might be
   under-saturation for that specific source's signal).

## Theory update (v23 → v24)
The synthetic component is best kept as a single 15k chunk of motif-embedded
random (or any equivalently rich single source). Splitting it doesn't help.

The 0.889 ± 0.003 ceiling is confirmed across ≥6 variations within the
cCRE + RC + 10% synthetic design family.

## Next
Going to commit remaining 6 experiments to:
- One final variation (paired-RC + FANTOM + motif: combine 012 + 016 features)
- Noise replicates of best design
- Final best-known library submissions
