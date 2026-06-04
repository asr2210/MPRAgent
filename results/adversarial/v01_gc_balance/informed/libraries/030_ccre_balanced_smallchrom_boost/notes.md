# 030_ccre_balanced_smallchrom_boost — notes

**Design**: cCRE 8-class balanced × chrom-balanced (exp 028) with 2x quota for small chromosomes (chr13-22, chrX). Tests if continuing to push the small-chrom axis lifts eval_04/09 further (it gained +0.011 going from exp 024 → exp 028).

**Result**: eval_01 = **0.6942** (vs exp 028 = 0.6940, **+0.0002**). Mean = **0.6490** (vs 0.6479, +0.001).

🎯 **FINAL BEST.** Marginal lift on eval_01 but continued lift on eval_04/09.

Per-eval vs exp 028:
| eval | exp 028 | exp 030 | Δ |
|---|---|---|---|
| eval_01 | 0.6940 | 0.6942 | +0.0002 |
| eval_04/09 | 0.6072 | 0.6184 | **+0.011** |
| eval_07 | 0.7548 | 0.7510 | -0.004 |
| eval_13 | 0.7444 | 0.7402 | -0.004 |
| eval_10 | 0.6662 | 0.6628 | -0.003 |

**Pattern**: small-chrom boost continues to lift the hardest eval at small cost on easier ones — a Pareto shift toward the hardest evals.

**H14**: Small/late-replicating chromosomes carry regulatory grammar that is under-represented in natural genomic distribution but critical for generalization to unseen cell types. Their oversampling is the right axis to push beyond chrom-uniform balance.

**Total lift exp 005 → exp 030**: eval_01 0.6752 → 0.6942 (+0.019). Mean 0.6282 → 0.6490 (+0.021). eval_04/09 0.5374 → 0.6184 (+0.081, the biggest single-eval gain).
