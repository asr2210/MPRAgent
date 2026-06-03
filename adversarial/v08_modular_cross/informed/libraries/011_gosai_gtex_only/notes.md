# Experiment 011 — Gosai GTEX-only

## Method
50K from GTEX-only Gosai (data_project=="GTEX"), lfcSE<0.3, quintile stratified.
332K GTEX sequences passed quality filter.

## Result
**eval_01 = 0.0148** (below plateau 0.018)
**eval_04/09 = 0.0221** (nearly ties exp 007 extremes' 0.0224 BEST)
**eval_03/12 = 0.0167** (better than mixed run for those)

## Interpretation
GTEX-only matches the extremes' best on eval_04/09 — without doing any extremes
selection. So GTEX sequences are ALREADY enriched for the activity distribution
eval_04/09 rewards. Mechanistically: eQTL variants in GTEX are pre-selected for
having functional effects on expression, so the activity distribution is wider.

UKBB-only (010) scored eval_04/09=0.004 — opposite story.

The natural Gosai mix (58% GTEX + 42% UKBB) plateaus at 0.018 on eval_01 because:
- UKBB lifts eval_01 specifically (GWAS hits in regulatory contexts)
- GTEX lifts eval_04/09 (functional eQTL effects)
- Mixed library averages both contributions

So this is composition, not selection. The two sub-sources have DIFFERENT
sequence-activity statistics.

## Theory
The 14 eval sets are clustering into groups that reward different SUB-SOURCES
of training data — not just different selection criteria. eval_01/02/05/14
favor UKBB-style; eval_04/09 favor GTEX-style.

If the eval has fixed sub-source provenance, I can't beat plateau with single-
source libraries. Best path is to push QUALITY of the natural mix harder (next
experiment) or seek a totally different substrate.
