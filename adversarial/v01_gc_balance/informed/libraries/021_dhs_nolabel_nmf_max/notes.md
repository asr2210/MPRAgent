# 021_dhs_nolabel_nmf_max — notes

**Design**: combine exp 005's no-label-topic filter with exp 020's NMF max-loading weighting.

**Result**: eval_01 = 0.6675 (vs exp 005 = 0.6752, -0.008). Notable lift on eval_04/09 (0.5579 vs 0.5374, +0.020).

**Conclusion**: combining the two doesn't beat uniform-within-filter for the primary metric. NMF weighting has a small but consistent eval_04/09 lift across exp 019/020/021 — likely indicating eval_04/09 specifically tests cell-type-specific regulatory grammar. But the cost is dropping eval_01 and others.

The NMF axis has informative substructure (different evals respond differently) but no global improvement in this pipeline.
