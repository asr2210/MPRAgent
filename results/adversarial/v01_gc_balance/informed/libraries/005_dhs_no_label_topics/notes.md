# 005_dhs_no_label_topics — notes

**Design**: 50k DHS uniformly sampled from the 12 NMF topics NOT aligned with labeling cells. Excluded: Myeloid/erythroid (K562), Digestive (HepG2), Cancer/epithelial (HepG2), Neural (SK-N-SH). Pool size: 2.61M of 3.59M (~73%).

**Prediction (H2)**: If model learns universal regulatory grammar, eval_01 stays near 0.66 (no significant drop). If labeling-cell-type alignment matters, eval_01 drops.

**Result**: eval_01 = **0.6752** (vs 0.6604 for full DHS uniform — actually *improved* by 0.015). Mean across 14 = 0.6311 (vs 0.6128 — improved 0.018).

**This is a strongly positive result for H2**. Removing data from the labeling cell types' aligned chromatin programs IMPROVED predictions on those very cells. The model trained on "other cell types' regulatory elements" predicts K562/HepG2/SK-N-SH activity *better* than the model trained on all elements including the matching ones.

**Why?** Two non-exclusive interpretations:
1. The 4 removed topics are large and partially redundant. Removing them gives the remaining 12 topics more relative weight, providing broader regulatory grammar coverage per training example.
2. Cell-type-specific topics may carry narrow regulatory signatures (e.g. lineage-specific TF motif arrangements) that overfit when present in training while not contributing universal grammar.

Either way, this directly answers the prompt's framing question: **yes, a library trained without test-cell-type-specific data can still predict on those cells**. The library *generalises*. This is the most theory-relevant result so far.

**Per-eval gains/losses**:
- eval_01,02,05,14: 0.660 → 0.675 (+0.015) — small uniform lift
- eval_03,12: 0.673 → 0.686 (+0.013)
- eval_04,09: 0.493 → 0.537 (+0.045) — biggest lift, on the previously-weak eval
- eval_07: 0.761 → 0.760 — flat
- eval_08: 0.116 → 0.121 — flat (still terrible)
- eval_10: 0.655 → 0.664
- eval_13: 0.752 → 0.751 — flat

The biggest lift is on eval_04/09 (which had moderate performance) — narrowed-down regulatory programs give *cleaner* training signal for moderate-activity prediction. eval_07/13 (motif-content-driven) plateau because they're already at the model's ceiling for natural regulatory sequences.

**Theory after 5 experiments (H2 → H2+)**:
> Library utility is driven by *diverse universal regulatory grammar*, not cell-type alignment. Removing data from labeling cell types' specific programs does not hurt; it can even slightly help by reducing overfitting to narrow lineage signatures. The library design strategy should optimise for *coverage of diverse natural regulatory contexts*, indifferent to whether those contexts are active in the labeling cell types.

**Next**: stack improvements. Use the "exclude label-aligned topics" winning pool, then add neutral random filler (which is known to be safe per `dhs_synth` baseline 0.7174 vs dhs_random 0.7089). Aim: ~0.70+.
