# 012_dhs_all_components — notes

## Design
3,125 DHS summit-centered 200bp windows from EACH of 16 NMF tissue
components (50k total). Opposite of 011 — maximum tissue-type breadth.

## Result
eval_01 = 0.0464  (vs dhs_specific=0.049, dhs_random=0.041, celltype_match=0.040)
mean_r = 0.0430

Beats dhs_random and celltype_match cleanly. Falls just short of dhs_specific.

## Interpretation
Tissue-breadth helps! Spreading evenly across 16 components beats focusing
on 4. But dhs_specific (weighted by 1/sqrt(numsamples)) is still better —
that weighting implicitly does breadth AND emphasizes rare sites.

This supports the theory: the eval rewards REPRESENTATION OF DIVERSE
REGULATORY DOMAINS. Specificity weighting and component balancing both
help. Cell-type-focused subsetting and signal-strength-focused subsetting
both hurt.
