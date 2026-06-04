# Experiment 003 — dhs80_synth20

## Design
40 000 uniform-random DHS + 10 000 i.i.d. uniform random synthetic.

## Hypothesis
eval_08 jumps (it likely contains synthetic sequences and only a model
trained on some non-genomic content can score it). eval_01 holds or
drops slightly.

## Result
| eval | exp 002 (100 % DHS) | exp 003 (80/20) | Δ |
|------|--------------------:|----------------:|------:|
| 01   | 0.5627 | 0.5273 | **−0.035** |
| 08   | 0.1663 | 0.2468 | **+0.081** |
| mean | 0.516  | 0.502  | −0.014 |

Direction matched (eval_08 up, eval_01 down) but net is **worse**.

## Interpretation
20 % synthetic is too much. The eval_01 cost (−0.035) exceeds the
eval_08 gain (+0.081 / 14 evals ≈ +0.006 averaged). And eval_08 only
moves from 0.17 to 0.25 — barely a third of the way to the
`random_uniform` baseline's 0.58 on eval_08.

Two readings:
1. **eval_08 is dominated by training distribution.** To score well on it
   you basically need a synthetic-majority library, which trashes every
   other eval. Pursuing eval_08 is a trap if eval_01 is the goal.
2. **DHS-derived sequences have a more efficient information-per-slot
   density than random.** Trading 10 k DHS slots for 10 k random
   sequences loses real predictive signal on the other 13 evals to chase
   one outlier eval.

## Theory update
- The "small OOD fraction helps generalization" hypothesis (item 4 from
  the initial theory) is partially wrong at the level I tested. 20 % is
  too much; maybe 1–5 % would be cost-neutral. But the bigger lesson is
  that eval_08 is in a different regime, and chasing it costs the rest.
- For maximizing eval_01, **dilute the DHS pool only if the substitute
  has higher predictive density than uniform random**. Pure random has
  near-zero predictive density.

## Numbers
mean_r averaged across 14 evals: 0.502
eval_01: 0.5273 (vs 0.5627 in exp 002)
eval_08: 0.2468 (vs 0.1663 in exp 002)
time_s: 30.5

## Next
Focus on the DHS portion. Try quality-filtering: drop noisy 1-biosample
peaks and very-low-signal calls (which dilute the regulatory signal).
Filtered pool size: ~1.18 M elements at (mean_signal > 0.5 &
numsamples >= 3) — plenty to sample 50 k from.
