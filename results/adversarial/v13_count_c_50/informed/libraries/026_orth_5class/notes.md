# Experiment 026 — drop CA-CTCF too (5-class cCRE)

## Result
eval_01 = 0.5747 (−0.0015 vs 025's 0.5762).

Dropping CA-CTCF hurt. The CTCF-bound class contributes meaningful
sequence information (likely CTCF motif itself, plus surrounding
insulator-context grammar). Pruning stops being beneficial after
dELS + CA.

## Numbers
mean_r: 0.555
eval_01: 0.5747
eval_07: 0.6083 (slight loss vs 025)
eval_08: 0.1717 (slight lift — more PLS+pELS per class)

## Conclusion
**The 6-class no-dELS-no-CA cCRE (025) is the sweet spot.** Further
class pruning over-trims and loses CTCF-related signal.

## Next
Stop pruning. Try complementary refinements:
- 027: seed stability check (rerun 025 with seed=1)
- 028: stricter numsamples cap (≤2 instead of ≤5)
- 029: small CpG slice on top of 025
