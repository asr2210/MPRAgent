# Experiment 019 — pure 50k orth-DHS

## Result
eval_01 = 0.5205 (vs 002 full DHS uniform: 0.5627; **−0.042**).

Pure orthogonal-DHS is WORSE than full DHS uniformly sampled.

| eval | 002 (full DHS) | 019 (pure orth DHS) | Δ |
|------|---------------:|--------------------:|--:|
| 01   | 0.5627 | 0.5205 | −0.042 |
| 04/09| 0.5405 | 0.4278 | −0.113 |
| 07   | 0.5974 | **0.6267** | +0.029 (best!) |
| 08   | 0.1663 | **−0.041** | −0.207 (negative!) |
| 13   | 0.5771 | 0.6099 | +0.033 |

## Interpretation
**The orth-DHS pool's value is COMPLEMENTARY to cCRE, not stand-alone.**

By itself, orth-DHS is missing the high-confidence regulatory regions
(PLS/pELS/active dELS) that drive eval_01 and especially eval_08.
eval_08 going NEGATIVE confirms it: orth-DHS has no CpG-rich content,
so a model trained on it actively *misranks* eval_08 sequences.

But orth-DHS gives huge eval_07 (+0.029) and eval_13 (+0.033) lift —
broad cell-type sequence diversity is exactly what those evals reward.

## Theory update
- **Library composition is a balancing act.** The "right" amount of
  each prior depends on which eval you optimize.
- **Orth-DHS without cCRE is BAD for eval_01**. The cCRE half provides
  high-confidence regulatory grammar essential for eval_01's task.
- **Optimal mix is ~25 k orth-DHS + 25 k cCRE class-bal (015)**:
  cCRE gives the regulatory grammar; orth-DHS adds non-cCRE diversity.

## Numbers
mean_r: 0.476 (low because eval_04/08 crashed)
eval_01: 0.5205
eval_07: 0.6267
eval_08: −0.041 (worst ever — confirms eval_08 wants CpG-rich)
time_s: 28

## Next
Don't go pure-orth. Stay with the mix. Try refining the orth-DHS
side: cell-type-specific subset (numsamples ≤ 5).
