# 007_dhs_nolabel_numsamples — notes

**Design**: Same pool as exp 005 (12 non-label-aligned NMF topics), but sample weighted by `numsamples` (number of biosamples each DHS is detected in).

Expected mean numsamples in sample after weighting: 254 (vs ~7 if uniform — so a 36x shift toward broadly accessible DHSs).

**Result**: eval_01 = 0.6185 (vs 0.6752 uniform from same pool). Mean across 14 = 0.5798 (vs 0.6311 uniform).

**Weighting hurts.** Across-the-board lower, biggest drops on eval_07/13.

**Why**: weighting by `numsamples` concentrates the library on broadly accessible / tissue-invariant DHSs. These are biologically robust but **homogeneous** — many similar housekeeping promoters / CTCF sites / ubiquitous TFBSs. The training set loses sequence/regulatory-context diversity. Cell-type-specific (low-numsamples) DHSs are the workhorse of learning because they carry varied regulatory programs.

**Theory update H3 → H4**:
> Within DHS, **diversity of regulatory programs** beats **biological robustness**. The model learns from contrast across many different regulatory contexts, not from repeated exposure to the same few "always-on" elements. Sampling that concentrates the library reduces diversity and hurts learning.

This is symmetric with the prior finding (exp 005): exclusion of dominant topics helped (more diversity); now weighting toward a narrow subclass (broadly-active) hurts (less diversity). The unifying principle is DIVERSITY.

**Next**: test if `mean_signal` weighting helps (different axis from `numsamples`). High `mean_signal` could correlate with strong activating elements specifically, providing cleaner activity labels per sequence. Predicts: small lift if signal-strength helps; flat if uniform already captures.
