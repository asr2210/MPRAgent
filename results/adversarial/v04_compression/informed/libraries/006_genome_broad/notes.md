# 006_genome_broad

## Design
50k random 200bp tiles from hg38 chr1–chr22 + chrX + chrY (all 24
chromosomes), weighted by chromosome length. seed=0.

Robustness check on 002 (which used only chr1/17/19/22).

## Result
- eval_01 = **0.4823** (vs 002's 0.4992; Δ -0.017)
- eval_07 = 0.5951 (vs 0.5985; Δ -0.003)
- eval_13 = 0.5905 (vs 0.6020; Δ -0.012)
- eval_04 = 0.4882 (vs 0.5226; Δ -0.035) — largest drop
- eval_08 = 0.0883 (vs 0.0916; Δ -0.003)
- mean over 14 ≈ 0.490 (vs 0.503; Δ -0.013)

## Interpretation
Broadening from 4 gene-rich chromosomes to all 24 chromosomes is
SLIGHTLY WORSE — by 0.013–0.017 r-points. Edge of noise but
consistent across most evals. Implies chr1/17/19/22 have a small
intrinsic advantage, plausibly:
- Higher gene density → more regulatory elements per random tile
- chr19 is the most gene-dense chromosome (especially GC-rich)
- chr17 also gene-rich

The advantage is small (<0.02). Probably not worth optimizing further
along this axis. My baseline is ~0.49–0.50 for "random real DNA" no
matter the chromosome subset. Lock that in.

## Implications
- Real-DNA random sampling has a ~0.50 plateau on eval_01 regardless
  of chromosome choice in this pipeline.
- The next experiment should test designs that could break this
  plateau: motif-implanted scaffolds, cCRE-enriched, activity-targeted,
  or cross-species real DNA.
- Don't expect chromosome-level optimization to add more than ~0.02.
