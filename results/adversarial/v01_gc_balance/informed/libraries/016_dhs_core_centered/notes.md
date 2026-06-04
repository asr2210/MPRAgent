# 016_dhs_core_centered — notes

**Design**: same filter as exp 005 (12 non-label-aligned topics), but window centered on `core_midpoint` instead of `summit`. Drops NaN-core rows.

**|core_midpoint - summit|** distribution: median 0, p90 23bp, max 1392bp. For most DHSs the two anchors coincide; even at the p90 the offset is well within the 200bp window.

**Result**: eval_01 = 0.6750 (vs exp 005 = 0.6752). **Flat — within noise.**

**Conclusion**: window anchor (summit vs core_midpoint) doesn't matter at this resolution. The 200bp window absorbs the small differences between the two coordinate definitions.
