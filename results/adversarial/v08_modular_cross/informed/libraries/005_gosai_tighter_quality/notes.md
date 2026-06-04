# Experiment 005 — Gosai tighter quality, per-cell stratification

## Method
- lfcSE < 0.3 filter (vs 0.5 in exp 004) → 599K pool
- 5x5x5 = 125 bins on (K562, HepG2, SKNSH) quintile activity, ~400/bin
- Top up with random fill

## Result
**eval_01 = 0.0180** (vs exp 004's 0.0181). Essentially same. Per-cell scores
more BALANCED but mean unchanged. Plateau confirmed.

## Implication
Tightening quality and refining stratification both flat. The bottleneck is
elsewhere — likely COMPOSITION of training set, not quality of labels.

Next: vary Gosai sub-source composition (UKBB vs GTEX vs CRE).
