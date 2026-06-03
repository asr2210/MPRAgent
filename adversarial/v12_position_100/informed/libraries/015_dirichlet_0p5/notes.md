# 015 — Dirichlet(0.5)

## What I built
50k sequences from Dirichlet(0.5). Alpha refinement test.

## Result
- eval_01 = 0.0768 (vs Dirichlet(0.3) = 0.0786, -0.0018)
- mean ≈ 0.0944 (vs 0.0954)
- eval_08 = 0.0726 (vs 0.0716, +0.0010) — slight improvement on noisiest eval
- All other evals slightly worse

## Interpretation
**Peak confirmed at alpha=0.3 for eval_01.** Full alpha sweep:
- 0.1 → 0.0752
- 0.3 → 0.0786 (PEAK)
- 0.5 → 0.0768
- 1.0 → ~0.0765 (exp 002 implicit)

Peak is sharp at alpha=0.3 for eval_01. Composition variance lever exhausted.

Interesting nuance: eval_08 likes alpha=0.5 better than 0.3 (different evals have
different alpha preferences). This suggests different evals weight composition
extremity differently — eval_01 favors more extreme, eval_08 favors slightly less.

## Hypothesis killed
H15: alpha > 0.3 is the peak → REJECTED. Peak is at 0.3.

## What to try next
Composition lever exhausted. Try strategy combination: 25k Dirichlet(0.3) + 25k shuffled
DHS. Shuffled DHS is the best biology variant (0.0754); maybe combining best synthetic
with best biology gets above 0.080.
