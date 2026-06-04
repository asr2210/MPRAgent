# Skill: Gosai 2024 MPRA dataset access

## Why it matters
v08's prepare.py uses an MPRA oracle (most likely Malinois) trained on the
Gosai 2024 K562/HepG2/SK-N-SH MPRA dataset. Training libraries must use
in-distribution Gosai sequences to get usable labels.

## Source data
- Paper: Gosai et al. 2024, Nature ("Machine-guided design of cell-type-targeting
  cis-regulatory elements")
- Boda2 codebase: https://github.com/sjgosai/boda2
- Table S2 MPRA dataset (798,064 rows):
  https://storage.googleapis.com/tewhey-public-data/CODA_resources/Table_S2__MPRA_dataset.txt

## Schema
Tab-separated, first row is header. Columns:
- IDs: variant identifier (e.g., "7:70038969:G:T:A:wC")
- chr: chromosome (1-22, X, Y as plain string, no "chr" prefix)
- data_project: source — UKBB (338k), GTEX (446k), CRE (14k)
- OL: oligo position
- class: phenotype label (e.g., "BMI,BFP" or "Depression_GP")
- K562_log2FC, HepG2_log2FC, SKNSH_log2FC: activity (log2 fold change)
- K562_lfcSE, HepG2_lfcSE, SKNSH_lfcSE: per-cell-type standard error
- sequence: 200bp ACGT (763,684 are exactly 200bp; rest 73-199bp — filter out)

## Splits (Gosai's own splits per Boda2)
- Train: chromosomes EXCEPT 7, 9, 13, 21, X
- Validation: 7, 13
- Test: 9, 21, X

In v08: chromosome-holdout filtering had no benefit (003 vs all-chrom 004),
suggesting v08 eval doesn't strictly respect Gosai chromosome split.

## Quality filtering insights (from experiments 003-004)
- 50K random Gosai (any chrom): eval_01 ≈ 0.0
- 50K Gosai with mean lfcSE < 0.5, activity-stratified: eval_01 ≈ 0.018
- Quality + activity stratification gives ~20x improvement
- Tighter lfcSE filtering should help further

## How to use
1. Read the file, skip rows where sequence length != 200 or non-ACGT chars
2. Compute mean_lfcSE = (K562_lfcSE + HepG2_lfcSE + SKNSH_lfcSE) / 3
3. Compute per-cell or mean activity for stratification
4. Filter by mean_lfcSE threshold (0.3-0.5)
5. Stratify-sample to ensure activity range coverage
