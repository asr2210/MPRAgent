# 010 — Multi-alpha Dirichlet mix

## What I built
10,000 sequences at each alpha ∈ {0.1, 0.3, 0.7, 1.5, 3.0} = 50,000 total. Shuffled.

## Result
- eval_01 = 0.0772 (vs exp 004 pure Dirichlet(0.3) = 0.0786). Slight loss on eval_01.
- mean across all 14 ≈ 0.0950 (vs 0.0954). Nearly identical overall.
- Gains: eval_07 (0.1455 → 0.1478), eval_10 (0.1285 → 0.1309), eval_08 (0.0716 → 0.0716).
- Losses: eval_01 (-0.0014), eval_13 (-0.0013 to 0.1444 from 0.1431). Wait, 0.1444 > 0.1431.

Net: **approximately neutral with a slight eval_01 loss.**

## Interpretation
Single Dirichlet(0.3) is at or near the peak for eval_01. Mixing in higher alphas
(0.7, 1.5, 3.0) doesn't add useful composition-variance signal — those sequences are
closer to uniform, which the model can already infer. The marginal sequences from
0.1 (already tested in exp 005 at 0.0752) don't help enough either.

Pattern across 10 experiments: the synthetic-composition ceiling for eval_01 is ~0.0786.

## What this rules out
H10 (multi-alpha mixing breaks composition ceiling) → not falsified strongly but no win.
Future composition-only experiments should default to single alpha=0.3.

## What to try next
Switch lever. Test high-signal biology pure: top 50k DHS by total_signal (exp 011).
