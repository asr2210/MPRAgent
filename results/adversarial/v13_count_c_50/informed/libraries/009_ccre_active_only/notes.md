# Experiment 009 — ccre_active_only

## Design
Uniform sample from cCRE classes with explicit active regulatory marks:
PLS (promoter-like), pELS (proximal enhancer-like), CA-H3K4me3. Pool
≈ 376 K.

## Hypothesis
Higher per-element regulatory info density → should beat broader cCRE
mixes. Tests "concentrate on active elements" intuition.

## Result — bad on eval_01, big on eval_08
| eval | 002 DHS uni | 006 cCRE CB | **009 active**|
|------|--------:|--------:|--------:|
| 01   | 0.5627 | 0.5637 | **0.4833** (−0.080)|
| 04/09| 0.5405 | 0.5841 | 0.5879 |
| 07   | 0.5974 | 0.5707 | **0.3900** (−0.181)|
| 08   | 0.1663 | 0.2333 | **0.3879** (+0.155)|
| 13   | 0.5771 | 0.5473 | **0.3668** (−0.181)|
| mean | 0.516  | 0.547  | 0.467  |

## Interpretation
**Confirmed**: shrinking to "highly active" elements destroys diversity-
dependent evals. eval_07 and eval_13 collapsed by ~0.18 each. eval_01
dropped by 0.08.

**Surprise**: eval_08 jumped to its highest level yet (0.388). Promoters
and proximal enhancers have *unusual base composition* (CpG-rich,
G-skewed) — far from average open chromatin. The model trained on these
sees more unusual sequence statistics, and apparently eval_08 contains
sequences with similarly unusual statistics.

This means **eval_08 is not really "synthetic-like"**. It's "atypical
genomic"-like. Promoters happen to look atypical compared to the average
DHS sequence because of CpG-island structure.

## Theory update
- **eval_08 = atypical genomic distribution.** Not random synthetic
  (which would give a flat sequence-statistic profile); not standard DHS
  (which gives a typical genome+regulatory profile). Promoters and
  similar low-complexity regulatory regions sit in the eval_08 wheelhouse.
  This is a *much* more interpretable model than "synthetic OOD."
- **Functional narrowness severely hurts eval_07 / eval_13.** These evals
  reward the model's ability to predict activity across a *wide spectrum*
  of regulatory contexts. Promoters alone don't teach it the spectrum.
- The "diversity dominates per-element density" rule is now reconfirmed
  twice (exp 004 and exp 009). The library should look like the
  regulatory grammar at large, not its highest-confidence subset.

## Numbers
mean_r: 0.467
eval_01: 0.4833
eval_08: 0.3879 (best so far)
time_s: 22.4

## Next
010 should test whether adding a small targeted PLS/pELS slice to the
008 mix can lift eval_08 without crashing the diversity-dependent evals.
e.g., 40 % DHS + 40 % cCRE class-balanced + 20 % PLS+pELS only.
