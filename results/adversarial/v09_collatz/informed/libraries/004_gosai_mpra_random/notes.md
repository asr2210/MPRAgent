# E004 — Gosai MPRA random sample

## Design
50k sequences sampled uniformly from the Gosai et al. 2024 MPRA dataset
(776K sequences measured in K562/HepG2/SK-N-SH at 200bp). This is almost
certainly the source dataset for `mpra_real` baseline (0.6026).

## Result
eval_01 = **0.3233** — barely above uniform DHS (E002=0.318), FAR below
the baseline `mpra_real` (0.6026).
Per cell: K562=0.149, HepG2=0.205, SKNSH=0.616.

## Critical interpretation
Gosai sequences (the source for `mpra_real`) score 0.32 in MY pipeline but
were reported as 0.60 in baseline. **Big gap**. Possible explanations:
1. Different prepare.py version — but I cannot verify
2. The baseline `mpra_real` may actually have used Gosai labels too (the
   "real" measurements), training a model with HIGH-quality labels, then
   evaluated on the same eval sets
3. My prepare.py "MPRA" may add more noise than the original Gosai labels

What I CAN observe:
- My pipeline produces consistent SKNSH ≈ 0.62, K562 ≈ 0.15, HepG2 ≈ 0.20
  across ALL my libraries (random, DHS, Gosai) — suggesting a structural
  ceiling per cell type
- K562 and HepG2 are systematically underlearned vs SKNSH
- The 14 eval sets cluster: identical scores for eval_01=eval_05=eval_14,
  eval_02=eval_06=eval_14, eval_04=eval_09, eval_11=eval_06, eval_03=eval_12

## Distribution of Gosai labels
- log2FC means: K562=0.51, HepG2=0.47, SKNSH=0.42 (similar)
- log2FC std: ~1.1-1.2 in each cell
- Median lfcSE: 0.19 (some sequences have very noisy labels)
- Cross-cell variance: mean 0.19, top 5% > 0.75

## Implications
Random sampling from Gosai is no better than random DHS. Need SMART selection:
- Sequences with strong, confident labels (high |log2FC|, low lfcSE)
- Maximize coverage of cell-type-specific activity patterns
- Avoid sequences with near-zero activity (noise-dominated labels)
