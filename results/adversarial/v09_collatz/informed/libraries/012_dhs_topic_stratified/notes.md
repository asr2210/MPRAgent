# E012 — DHS topic-stratified

50K = 3125 sequences per NMF component (16 components from DHS Index).

## Result
eval_01 = 0.3181. Per-cell: K562=0.144, HepG2=0.195, SKNSH=0.614.

## Interpretation
Essentially identical to E002 (DHS uniform random = 0.3179). Forcing equal
representation across the 16 chromatin programs gives no measurable benefit
over uniform sampling. Matches the baseline Table 1 finding
(dhs_stratified 0.7055 vs dhs_random 0.7089 — also indistinguishable).

Topic stratification is not a useful lever within DHS.
