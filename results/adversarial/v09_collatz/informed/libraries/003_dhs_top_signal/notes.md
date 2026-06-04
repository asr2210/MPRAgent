# E003 — DHS top mean_signal × log(numsamples)

## Design
Top 50k DHS sites by mean_signal × log(1+numsamples). Hypothesis was that
high-signal + multi-tissue sites would be the most measurably active under
real labels → higher SNR → better training.

## Result
eval_01 = **0.2661** — WORSE than uniform DHS-random (E002=0.318) and barely better than pure random (E001=0.244).
Per cell type: K562=0.141, HepG2=0.034, SKNSH=0.623.
HepG2 dropped from 0.19 (uniform DHS) to 0.03 — HUGE drop.

## Interpretation
Top-signal DHS sites are dominated by **constitutively active promoters**
(tissue-invariant, broadly accessible). These sites are:
- All ~similar high-activity → low variance for model to learn from
- Narrowed regulatory grammar → fewer distinct TF combinations
- Especially poor for HepG2 prediction (anti-selected against cell-type specific signal)

**Diversity beats SNR-maximization**. This contradicts my E002→E003 hypothesis.

## Theory update
- Pure SNR-maximization doesn't help; it narrows regulatory grammar
- The model needs sequences with VARIED activity, not all-high activity
- A library that all looks the same (high signal, similar context) = uninformative
