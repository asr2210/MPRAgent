# Experiment 005 — dhs_component_stratified

## Design
3 125 DHS elements per NMF component (16 × 3 125 = 50 000), uniform
within each. Forces balance across heavily imbalanced component pool
(Primitive/embryonic 626 K vs Stromal A 56 K).

## Hypothesis
After exp 004, the working theory said diversity dominates per-element
quality. Stratifying should improve eval_01 by upweighting under-
represented chromatin programs.

## Result
| eval | exp 002 (uniform) | exp 005 (stratified) |
|------|------------------:|---------------------:|
| 01   | 0.5627 | 0.5598 (−0.003)|
| 08   | 0.1663 | 0.1309 (−0.035)|
| 13   | 0.5771 | 0.5889 (+0.012)|
| mean | 0.516  | 0.513 |

Essentially a tie. eval_13 (which favors SEI in instructions baselines)
improved a touch, but eval_08 dropped and eval_01 was flat.

## Interpretation
**At 50k size, uniform DHS already captures enough cross-program
diversity.** Even Stromal A's 56K pool contributes ~770 sequences under
uniform sampling — enough for the model to see the program. Force-
balancing to 3 125/class doesn't move the needle.

This is a **negative result for the "diversity dominates" hypothesis at
the topic level.** Topic balance ≠ what matters. The real diversity axis
might be something more granular (motif identity, GC tier, etc.) that
NMF components only weakly correlate with.

## Refined theory
- Per-program (NMF topic) balance is approximately a free parameter at
  50k scale — pool imbalance is large enough that even the smallest
  program is well-represented.
- The next test is whether a different axis of diversity (functional
  regulatory class, not co-accessibility topic) actually matters.

## Numbers
mean_r averaged across 14 evals: 0.513
eval_01: 0.5598
eval_08: 0.1309
time_s: 12.4
