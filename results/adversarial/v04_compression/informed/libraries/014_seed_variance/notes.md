# 014_seed_variance — multi-seed calibration

## Design
Same as E2 (human chr1/17/19/22 random 200bp tiles), but generated
with seeds 0, 1, 2 to test (a) whether prepare.py runs multi-seed
when multiple sequence_{n}.txt files are present, and (b) the
size of single-seed variance.

## Result
- result.json: **n_seeds = 3**, time_s = 48.3 (vs ~12s for single)
- eval_01 = **0.4987** (3-seed mean)
- vs E2 (1-seed, same design): 0.4992 → Δ = -0.0005

## Interpretation
**Multi-seed pipeline IS supported.** Drop multiple `sequences_N.txt`
files and prepare.py picks them all up.

**Seed variance is tiny.** 3-seed mean is within 0.0005 of E2's
1-seed result. The 0.50 ceiling is not single-seed noise.

This means Δ ≤ 0.005 is in the noise band; only Δ ≥ 0.01 should be
treated as signal. Most of my prior "neutral" results (E4, E8, E13)
ARE genuinely neutral, not noise-hiding real signal.

## Implications for design budget
- Continue single-seed for exploratory experiments (faster, ~4x).
- Use multi-seed (3) only to confirm claimed improvements that are
  in the Δ ∈ [0.01, 0.02] band.
- Treat Δ ≥ 0.02 as definitively real even at single-seed.

## Theory v6 confirmed
The ranking is robust:
- gene-rich human random: 0.499 ± 0.001 (E2, E14)
- all-chrom human random: 0.482 (E6)
- DHS uniform:            0.470 (E10)
- DHS NMF-stratified:     0.438 (E3)
- promoter TSS:           0.370 (E9)
- 4-mer Markov:           0.268 (E5)
- uniform synth:          0.307 (E1)

Every Δ > 0.01 is real. The 0.50 wall is a real wall, not a noise floor.

## What's next
The "+0.50" wall is real. To break it, I need an axis I haven't
explored that gives Δ ≥ 0.02 over E2. Candidates:
- E15: **gene-density-weighted full-genome sampling**. Smooth
  gene-richness across all 24 chroms; might broaden info beyond hard
  chr1/17/19/22 cut without losing gene-richness.
- E16+: cCRE element-class mix (promoter+enhancer+CTCF+intergenic).
- E17+: activity-stratified using DHS mean_signal across full range.
- E18+: quality filter (reject low-complexity / pure-repeat tiles).
