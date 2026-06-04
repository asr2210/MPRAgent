# Experiment 024 — GTEX lfcSE<0.7 random, seed=123 (variance check)

## Result
**eval_01 = 0.0157** vs seed=42's 0.0222. Δ=0.0065.

## Interpretation — CRITICAL
Huge seed variance. The "breakthrough" at 0.0222 is partly luck of the draw.
True expected eval_01 for GTEX lfcSE<0.7 random is ~0.018-0.019, with the 0.022
being a high-side outcome.

## Implication
Random sampling has ~0.005 std at 50K from ~400K. Most of my fine-grained
"differences" between exp 017-023 are within this variance.

**Robust GTEX advantage:** all GTEX configs (017,019,020,022,023,024) score
0.014-0.022 vs mixed Gosai's ~0.017 ceiling. The GTEX edge is real but small.

## Path forward
To beat sampling noise, need DETERMINISTIC selection criteria rather than
random. Try:
- Take all GTEX with |mean activity| > threshold (effect-magnitude filter)
- Take top 50K GTEX by composite quality+signal score
- Or accept current best as the answer and document
