# Experiment 021 — orth-DHS + cCRE + CpG-rich supplement

## Design
22 k orth-DHS (numsamples ≤ 5) + 22 k cCRE class-bal +
6 k PLS+pELS supplement (3k each).

Tests whether a CpG-rich slice recovers eval_08 without crashing
eval_01.

## Result — eval_08 recovers, eval_01 holds
| eval | 015    | 020    | **021** | Δ vs 020 |
|------|-------:|-------:|--------:|---------:|
| 01   | 0.5736 | 0.5745 | **0.5746** | +0.0001 |
| 04/09| 0.5535 | 0.5509 | **0.5711** | **+0.020** |
| 07   | 0.6131 | 0.6168 | 0.6026 | −0.014 |
| 08   | 0.1492 | 0.1404 | **0.1836** | **+0.044** |
| 13   | 0.5924 | 0.5968 | 0.5813 | −0.016 |
| mean (14 evals) | 0.560 | 0.555 | **0.558** | +0.003 |

## Interpretation
**Best of both worlds (almost)**. Adding 6 k PLS+pELS:
- ties eval_01 with 020 (0.5746)
- recovers eval_08 substantially (+0.044)
- recovers eval_04/09 (+0.020)
- mildly loses eval_07/13 (−0.014, −0.016)
- mean_r lifts by +0.003 overall

**Why this works**: Orth-DHS provides cell-type-specific diversity
(lifts HepG2/SK-N-SH). cCRE class-bal provides regulatory grammar
breadth. PLS+pELS provides CpG-rich content (eval_08-relevant).
The 22/22/6 split spends most of the budget on the 015 winning
strategy while reserving 12 % for eval_08 rescue.

eval_07 and eval_13 dip because they reward broad cell-type
diversity, and the 6 k slice toward CpG-rich CC steals from that
diversity budget. Trade-off is mild though.

## Theory update
- **Mean_r can be improved beyond what eval_01 alone shows**.
  021 has same eval_01 as 020 but higher mean_r. If the harness's
  final score is a mean across evals, this is real progress.
- **Eval_08's atypical-CpG hypothesis is confirmed**. A 12 % CpG-
  rich slice gives a 0.044 lift on eval_08. The relationship is
  roughly linear (009 was 100 % CpG-rich → 0.388; 021 is 12 % → 0.183).

## Numbers
mean_r: 0.558
eval_01: 0.5746 (tied with 020)
eval_04/09: 0.5711
eval_08: 0.1836
eval_13: 0.5813

## Next
Try pushing CpG slice larger (8-10 k) to test where the eval_01
plateau breaks. Also explore cCRE class rebalance: reduce dELS
(huge but generic distal enhancers), boost rare classes (TF, CA-TF).
