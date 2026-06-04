# 024_ccre_class_balanced — notes

**Design**: Equal samples per of 8 cCRE classes (6250 each). Window centered on midpoint, ±100bp. Natural class proportions (exp 022) gave dELS 62.5% of draws; this forces 12.5% each.

**Result**: eval_01 = **0.6921** (vs exp 022 = 0.6827, **+0.0094**). Mean = 0.6460 (vs 0.6381, **+0.008**).

🎯 **NEW BEST** — second consecutive lift, larger than the cCRE-vs-DHS jump.

Per-eval lift over exp 022:
| eval | exp 022 | exp 024 | Δ |
|---|---|---|---|
| eval_01 | 0.6827 | 0.6921 | +0.009 |
| eval_02 | 0.6835 | 0.6930 | +0.010 |
| eval_03 | 0.6915 | 0.6995 | +0.008 |
| eval_04/09 | 0.5704 | 0.5966 | **+0.026** |
| eval_07 | 0.7554 | 0.7576 | +0.002 |
| eval_10 | 0.6660 | 0.6671 | +0.001 |
| eval_13 | 0.7451 | 0.7479 | +0.003 |

**Every eval improves.** The biggest gain is on eval_04/09 (the hardest evals), suggesting the rare classes (PLS=2%, CA-TF=1%, CA-H3K4me3=3%, TF=4%, CA-CTCF=5%) carry distinctive regulatory grammar that natural-proportion sampling under-samples.

**H11**: Class balancing across regulatory program types > natural proportions. The dELS-dominated natural mix concentrates too much variance on enhancer-like grammar; promoters (PLS/pELS) and CTCF/TF anchors contribute orthogonal signal. Diversity should be measured by *coverage of distinct mechanisms*, not by *frequency in the genome*.

**Next**:
- Exp 025: cCRE excluding PLS only (test if PLS specifically helps or if it's the balanced-rare-classes effect more broadly)
- Exp 026: union of cCRE class-balanced + DHS no-label (cross-source diversification on top of class balance)
- Exp 027: cCRE class-balanced with stricter per-class oversampling of rarest classes
