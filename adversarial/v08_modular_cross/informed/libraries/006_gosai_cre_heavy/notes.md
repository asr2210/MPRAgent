# Experiment 006 — CRE-heavy blend

## Method
- All 14K Gosai CRE-class (no SE filter — keep all)
- 36K UKBB+GTEX with lfcSE<0.3, per-cell quintile stratified

## Result
**eval_01 = 0.0105.** WORSE than the 0.018 plateau.

## Interpretation
CRE-heavy hurts. Either:
- CRE has poorer measurement quality (no SE filter applied) → noisy labels
- The eval distribution is variant-centric (UKBB/GTEX-like), not regulatory-
  element-centric, so spending 14K slots on CRE wastes capacity
- Both

The drop of ~0.007 from baseline tells me variant-centric Gosai is the better
substrate. CRE-class is the WRONG composition direction.
