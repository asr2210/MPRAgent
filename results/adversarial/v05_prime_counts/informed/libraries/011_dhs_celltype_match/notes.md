# 011_dhs_celltype_match — notes

## Design
50k DHS summit-centered 200bp windows, 12,500 from each of 4 components that
match the three measured cell types:
- Myeloid / erythroid → K562 (erythroleukemic)
- Cancer / epithelial → HepG2 (HCC)
- Digestive → HepG2-relevant
- Neural → SK-N-SH (neuroblastoma)

## Hypothesis
If the eval rewards activity in K562/HepG2/SK-N-SH, then DHS sites from
those exact components should give cleanest/strongest signal in MPRA and
beat agnostic DHS sampling.

## Result
eval_01 = 0.0401 (WORSE than dhs_random 0.041, dhs_specific 0.049)
mean_r = 0.0379
eval_08 = 0.0454 (modestly above bio baseline but below random/composite)

K562_r=0.034, HepG2_r=0.038, SKNSH_r=0.049 — interestingly SKNSH improved
slightly but K562 and HepG2 are below dhs_specific levels.

## Interpretation
Cell-type matching HURTS. The eval is NOT preferentially scoring sequences
that "look like" K562/HepG2/SK-N-SH regulatory DNA. It actually rewards
broader regulatory diversity. This rules out a "match the cell types"
strategy entirely. Specificity (1/sqrt(numsamples)) remains the best DHS
weighting — it diversifies across all 16 tissues, not just three.

## Next probe
Try k-mer-balanced synthetic or repeat-pattern probes. The cell-type
direction is dead.
