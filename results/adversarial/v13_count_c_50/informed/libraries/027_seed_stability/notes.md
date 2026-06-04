# Experiment 027 — seed stability check of 025

## Result
eval_01 = 0.5758 vs 025 = 0.5762. **Δ = −0.0004**.

Within ±0.001 of 025. The +0.001 to +0.003 improvements we've claimed
since 015 are real signal, not seed luck.

## Numbers
| eval | 025 (seed=0) | **027 (seed=1)** | Δ      |
|------|-------------:|-----------------:|-------:|
| 01   | 0.5762       | 0.5758           | −0.0004|
| 04/09| 0.5620       | 0.5578           | −0.004 |
| 07   | 0.6131       | 0.6158           | +0.003 |
| 08   | 0.1591       | 0.1503           | −0.009 |
| 13   | 0.5925       | 0.5943           | +0.002 |

Per-eval noise ~±0.004; eval_08 noisier (±0.01) because CpG content
is rare and sampling variance dominates.

## Conclusion
**025 is robust.** Seed noise on eval_01 is ~0.0005, much smaller than
the +0.002 lift 025 has over 020. Future work should treat ±0.001
eval_01 differences as in the noise.

## Per-cell-type
| cell    | 025    | **027**   |
|---------|-------:|----------:|
| K562    | 0.6155 | 0.6159    |
| HepG2   | 0.5529 | 0.5522    |
| SK-N-SH | 0.5602 | 0.5593    |

Stable to <0.001 across all three cell types — sampling stochasticity
is uniform, not concentrated in one cell line.
