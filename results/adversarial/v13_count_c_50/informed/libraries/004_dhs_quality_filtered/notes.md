# Experiment 004 — dhs_quality_filtered

## Design
Uniform random over DHS elements with `mean_signal > 0.5` AND
`numsamples >= 3` (1.18 M / 3.59 M survive).

## Hypothesis
Filtering noise should improve eval_01. Single-biosample peaks are
noisy, so removing them gives cleaner training signal.

## Result
| eval        | exp 002 (all) | exp 004 (filtered) | Δ |
|-------------|--------------:|-------------------:|------:|
| 01          | 0.5627 | 0.5127 | **−0.050** |
| 08          | 0.1663 | 0.2765 | +0.110 |
| mean        | 0.516  | 0.491  | −0.025 |
| K562 (01)   | 0.6073 | 0.6105 | +0.003 |
| HepG2 (01)  | 0.5273 | 0.4645 | **−0.063** |
| SKNSH (01)  | 0.5536 | 0.4632 | **−0.090** |

Filtering **hurt** eval_01 by 5 points. Worse, the regression was almost
entirely in HepG2 and SK-N-SH — K562 was unchanged.

## Interpretation
This is a clean diagnosis. "Quality" filtering by signal strength
preferentially keeps tissue-invariant / high-coverage elements (those
hit the signal threshold by being broadly accessible). That tilts the
training set toward K562-like cells, where the model has the most
collateral data. HepG2 and SK-N-SH suffer because cell-type-specific
DHS for those tissues tend to have lower signal / fewer biosamples and
were filtered out.

The "noise" I was trying to remove was actually carrying real cell-type-
specific information for the underrepresented cell types.

## Theory update — major

**At 50k library size, diversity dominates per-element quality.** A noisy
DHS call in a rare cell type contains information the model can use; a
high-confidence call in a tissue-invariant region is redundant with
many other elements already in the library.

This is the *opposite* of what I expected after exp 001. The reconciled
view: 1-biosample peaks weighted heavily (exp 001) are too noisy, but
filtering them out *entirely* (exp 004) costs too much diversity. The
sweet spot is uniform sampling over the full pool — let the law of
large numbers handle the noise.

This also predicts that for generalization to *unseen* cell types,
diversity of regulatory programs matters more than confidence in any
single element. The library should look like the chromatin universe, not
the cleanest 30 % of it.

## Numbers
mean_r averaged across 14 evals: 0.491
eval_01: 0.5127
eval_08: 0.2765
time_s: 12.4

## Next
Push diversity further: equal-per-component stratified sampling
(n/16 per NMF topic). The pool is very unbalanced (Primitive/embryonic
has 626 K elements, Stromal A only 56 K), and uniform sampling
inherits that imbalance. Force-balancing should oversample the small
components — directly testing the "diversity dominates" hypothesis.

If diversity is really king, equal-stratified should beat uniform
(reverse the strategies.md Table 1 ordering where `dhs_topic` >
`dhs_stratified`). If uniform still wins, the prior NMF imbalance
actually reflects information density and shouldn't be flattened.
