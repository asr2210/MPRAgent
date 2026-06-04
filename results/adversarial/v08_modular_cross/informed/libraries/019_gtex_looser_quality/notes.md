# Experiment 019 — GTEX looser quality (lfcSE<0.7) random — NEW BEST

## Method
GTEX-only with lfcSE<0.7 (409K — almost all of GTEX's 446K passes). Random 50K.

## Result
**eval_01 = 0.0222** — best yet. +17% over 017 (0.0190).
**eval_02 = 0.0225** best ever.
**eval_04/09 = 0.0222** ties extremes' best.
**Mean across 14 evals = 0.0164** — best aggregate.

## Interpretation
DIVERSITY beats QUALITY. The lfcSE<0.5 → <0.7 step allowed in many more
sequences (332K → 409K, adds 77K candidates) with mildly noisier labels.
Net effect: the model sees more sequence-activity examples, even if some are
noisier. Sequence-space coverage matters more than per-label cleanliness.

This inverts the earlier conclusion. Within GTEX, the optimal is:
- Loose quality filter (lfcSE<0.7) — barely filtering anything
- Random sampling (no stratification)
- Sub-source selection (GTEX-only)

## Theory update — the THREE KNOBS
1. **Sub-source**: GTEX-only > Gosai-mixed > UKBB > CRE
2. **Quality**: a U-curve. Too noisy fails. Too clean removes informative
   variation. Sweet spot ~lfcSE<0.7 for GTEX.
3. **Sampling shape**: random/natural > stratified > extremes/selection

## Next direction
- GTEX with NO filter (random 50K from 446K) — push diversity to max
- Test if lfcSE<0.7 helps mixed Gosai too (parallel test)
- Try GTEX + small CRE/UKBB augmentation
