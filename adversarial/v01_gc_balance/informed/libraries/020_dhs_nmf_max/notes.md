# 020_dhs_nmf_max — notes

**Design**: same as exp 019 but weight = per-DHS MAX NMF loading (instead of sum). Upweights DHSs strongly aligned with one specific topic.

**Result**: eval_01 = 0.6641 (vs exp 005 = 0.6752, -0.011; vs exp 019 sum = 0.6481, +0.016).

Max interpretation slightly better than sum, but both worse than uniform.

**Conclusion**: in this pipeline, neither NMF interpretation closes the gap to the baseline `dhs_topic` (0.7232). The information in NMF loadings doesn't translate to a measurable lift in MY prepare.py. Most likely a different surrogate / scoring model.

**Notable per-eval**: NMF weighting slightly LIFTS eval_04/09 (0.5409 vs 0.5374 uniform) — the hardest eval. Small effect, but consistent across exp 019/020/021. Suggests cell-specificity weighting helps a specific eval set, just hurts most others.
