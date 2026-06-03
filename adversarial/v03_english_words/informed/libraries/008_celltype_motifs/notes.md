# 008 — Cell-type-targeted JASPAR motifs

## Plan
Restrict the JASPAR motif pool to cell-type-relevant TFs (K562 erythroid,
HepG2 hepatocyte, SK-N-SH neural) plus a handful of universal regulators
(SP1, CTCF, NF-Y, AP-1, ETS, etc). Insert 3 per sequence in random uniform
backbone. Match exp 006's design but with focused motif pool (289 motifs vs 2344).

## Result
**eval_01 = 0.4283.** New best. Beats exp 006 (0.4252) by +0.003.
K562: 0.596 (↑ from 0.591), HepG2: 0.629 (↑ from 0.622), SK-N-SH: 0.060 (≈).

## What this teaches
- Cell-type-targeted motifs do help K562 and HepG2. The surrogate appears to
  reward sequences enriched for the right TFs.
- SK-N-SH is stuck at 0.06 even with neural-motif enrichment. Strong evidence
  that **SK-N-SH performance has a structural floor at ~0.06 in this pipeline,
  not addressable via library design alone**. The SK-N-SH oracle might be a
  noisy or weak target.
- Improvement is small (+0.003), suggesting limited room for further gain from
  motif-only design.

## Theory T7
- Cell-type-relevant motifs are the lever for K562/HepG2 (each adds ~1% per
  axis of targeting).
- SK-N-SH is a constant penalty; ignore it as an optimization target and focus
  on K562/HepG2 + mean_r.
- The MY-pipeline ceiling appears to be ~0.43-0.45 for mean_r. The path to
  higher numbers may not exist with simple motif insertion.

## Next
Test motif DENSITY with cell-type-only pool: 5 motifs/sequence (vs 3).
If higher density helps with focused pool: exp 009 > 0.43.
If saturation: exp 009 ≈ 0.42-0.43.

If saturated, pivot to mixture libraries (motif + DHS + random) or test
whether MORE training data per category (e.g., a 25k pure K562 motif library
+ 25k pure HepG2 motif library) outperforms uniform mixing of all three.
