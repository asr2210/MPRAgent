# 023_ccre_intersect_dhs_nolabel — notes

**Design**: cCRE all-classes filtered to those whose midpoint falls within a no-label-aligned DHS interval. 1.47M of 2.35M cCREs retained (62.6%). Uniform 50k sample.

**Result**: eval_01 = 0.6796 (vs exp 022 = 0.6827, **-0.003**). Mean = 0.6330 (vs 0.6381, -0.005). Slight decline.

**Conclusion**: combining topic exclusion (via DHS overlap) with cCRE selection slightly HURTS — removes 880k cCREs that were apparently useful. The 4 label-aligned NMF topics map to DHSs that overlap PRODUCTIVE cCREs; removing them costs more than the label-confound saves.

**H10 refinement**: cCRE's curation already integrates a more informative signal than DHS-NMF-topic exclusion. Adding the DHS topic filter on top is redundant or harmful.

**Best library remains exp 022 (cCRE all-classes uniform, 0.6827).**

**Next**: ablate cCRE classes directly to see which contribute most (vs the topic-overlap filter which is indirect).
