# 015_dhs_inv_numsamples — notes

**Design**: DHS no-label-aligned pool, sampled with weight ∝ 1/(1+numsamples). Favors cell-type-specific DHSs (detected in few biosamples). Expected mean numsamples after inverse-weight = 2.82 vs uniform mean ~26.14.

**Motivation**: re-reading baseline recipes — `dhs_topic` (0.7232 best baseline) samples proportional to NMF topic loadings, which UPweights cell-type-specific DHSs. Exp 007's `numsamples`-weighted (favoring broad DHSs) failed (0.6185). Predicted inverse-weight should help.

**Result**: eval_01 = 0.6598 (vs exp 005 uniform = 0.6752, exp 007 numsamples-weighted = 0.6185). Mean = 0.6097.

**Surprise**: inverse weighting also hurts, just less than forward weighting. The numsamples axis is asymmetric — favoring broad hurts MORE than favoring specific — but uniform still wins. The baseline `dhs_topic` lift presumably comes from full NMF topic loadings (continuous per-topic weights), not just the inverse numsamples proxy.

**Theory update H6 → H7**:
> On the numsamples (cell-specificity) axis, uniform sampling within a filtered pool still wins. The `dhs_topic` baseline's lift comes from a richer signal (per-DHS NMF topic loadings) not available in the simple `component` dominant-topic label. Without that signal, I cannot easily replicate the dhs_topic recipe.

**Implication**: stop chasing the 0.7232 baseline gap via re-weighting. The instrument I'd need (per-DHS topic loading matrix) isn't in the file I have. Move to orthogonal axes: window placement (core_midpoint vs summit), multi-view per source, peak-width filter.
