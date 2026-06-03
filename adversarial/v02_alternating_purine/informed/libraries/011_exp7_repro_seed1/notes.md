# 011_exp7_repro_seed1 — reproducibility test

**Composition:** identical to exp 007. Only change: SEED=1 throughout.

**Result:** mean_r(14) = **0.1498**, eval_01 = 0.1521. vs exp 007 (0.1532, 0.1592): -0.003 / -0.007.

**Per-eval delta (exp 007 → exp 011, both same recipe):**
| eval | exp 007 (seed 0) | exp 011 (seed 1) | Δ |
|---|---|---|---|
| eval_01 | 0.1592 | 0.1521 | -0.007 |
| eval_06/11 | **0.2283** | **0.1603** | **-0.068** |
| eval_07 | 0.1311 | **0.1863** | **+0.055** |
| eval_04/09 | 0.1357 | 0.1647 | +0.029 |
| eval_13 | 0.1266 | 0.1396 | +0.013 |
| mean14 | 0.1532 | 0.1498 | -0.003 |

**Findings:**
1. **mean14 is stable across seeds** at ~0.15 (Δ = 0.003). The aggregate is robust.
2. **Per-eval r is HIGHLY seed-sensitive** — eval_06/11 swung 0.068, eval_07 swung 0.055.
3. Exp 007's eval_06/11 = 0.228 was a seed-level lucky shot. The recipe's true eval_06/11 ≈ 0.190 (avg).
4. Exp 007's eval_07 = 0.131 was unlucky. True eval_07 ≈ 0.16.

**Implications for optimization:**
- I cannot trust per-eval gains <0.05 — those are within seed noise.
- mean_r(14) gains <0.005 are within noise; I need interventions with predicted lift ≥0.01.
- The "true" multi-source recipe baseline is ~0.151 mean_r.
- Need qualitatively different interventions, not subtle reweightings.

**Top reproduced result:** Multi-source recipe (exp 007/011 family) at mean_r ~0.151. Hard to beat without a structural change.
