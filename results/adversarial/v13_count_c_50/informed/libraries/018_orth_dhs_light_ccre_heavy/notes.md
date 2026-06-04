# Experiment 018 — orth-DHS light + cCRE heavy (15/35)

## Result
eval_01 = 0.5734 (tie with 015's 0.5736). The other-side ratio gives
the same answer. The orth-DHS+cCRE mix is robust in the 15-25k orth
range.

| eval | 015 (25/25) | **018 (15/35)** | 016 (35/15) |
|------|------------:|----------------:|------------:|
| 01   | 0.5736 | 0.5734 | 0.5689 |
| 04/09| 0.5535 | 0.5672 | 0.5279 |
| 07   | 0.6131 | 0.6036 | 0.6250 |
| 08   | 0.1492 | 0.1787 | 0.0972 |
| 13   | 0.5924 | 0.5814 | 0.6058 |

## Interpretation
015 and 018 statistical tie on eval_01. The ratio is flexible
between 15/35 and 25/25 (cCRE-heavy side). Beyond 25k orth-DHS, eval_01 falls.

Looking at the K562/HepG2/SK-N-SH split on eval_01:
- 015: 0.6131 / 0.5485 / 0.5592 = 0.5736
- 018: 0.6194 / 0.5456 / 0.5551 = 0.5734

018 has slightly higher K562 (0.6194 vs 0.6131) — more cCRE rescues
K562. But it gives back some HepG2/SK-N-SH lift. Net wash on eval_01.

eval_04/09 lifts to 0.5672 in 018 — closer to 012's 0.5772. More cCRE
clearly helps eval_04/09.

## Numbers
mean_r: 0.557
eval_01: 0.5734
eval_04: 0.5672
eval_08: 0.1787 (less bad than 015's 0.1492)

## Next
018 doesn't beat 015 on eval_01, so 25/25 stays the orth-DHS+cCRE
optimum. Try a more targeted orth-DHS variant — cell-type-specific
subset (low numsamples).
