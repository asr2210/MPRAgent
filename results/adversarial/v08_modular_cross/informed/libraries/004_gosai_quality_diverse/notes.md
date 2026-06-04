# Experiment 004 — Gosai quality + activity-stratified

## Method
- Filter Gosai 798K → 697K with mean lfcSE < 0.5
- Stratify by mean activity across cell types into 5 quintile bins
- Sample 10K per bin → 50K total

## Result
**eval_01 = 0.0181** (up from 0.0000 in exp 003 with random Gosai).
Most eval sets now in 0.01-0.02 range. First clear positive signal.

## Interpretation
Quality filtering and activity stratification provide ~20x improvement over
random Gosai. Two mechanisms likely at play:
1. **Cleaner labels**: lower-SE sequences give the prepare.py oracle (presumed
   Malinois-like) more reliable inputs, so per-example labels are less noisy.
2. **Activity diversity**: stratifying ensures full dynamic range coverage, so
   the model sees the full "active" and "inactive" regimes.

## Theory update
The prepare.py oracle is likely an MPRA-trained model (almost certainly trained
on Gosai). It gives accurate labels for in-distribution Gosai sequences and
arbitrary labels for OOD sequences (DHS, cCRE, random). The trained surrogate
model needs:
- Sequences the oracle labels precisely
- Full activity range coverage
- Diverse sequence content

For experiment 005: tighter quality (lfcSE < 0.3) + per-cell-type
stratification (each cell type has its own dynamic range, useful for learning
cell-specific activity).
