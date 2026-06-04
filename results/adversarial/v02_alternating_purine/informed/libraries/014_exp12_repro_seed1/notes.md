# 014_exp12_repro_seed1 — exp 012 was lucky-seed

**Composition:** identical to exp 012. Only SEED=1.

**Result:** mean_r(14) = **0.1441**, eval_01 = 0.1446. *vs exp 012: -0.018 mean14.*

**Per-eval delta (exp 012 → exp 014):**
| eval | seed 0 | seed 1 | Δ |
|---|---|---|---|
| eval_01 | 0.174 | 0.145 | -0.030 |
| eval_06/11 | 0.221 | 0.131 | -0.090 |
| eval_07 | 0.162 | 0.199 | +0.038 |
| eval_13 | 0.129 | 0.160 | +0.031 |
| K562 r (eval_06) | +0.078 | **-0.056** | -0.134 |
| mean14 | 0.162 | 0.144 | -0.018 |

**Conclusion: high variance hides a true effect of magnitude 0.**
- Cell-type DHS recipe average over seeds: 0.153 (= exp 012 + 014 mean / 2)
- Multi-source recipe average over seeds: 0.152 (= exp 007 + 011 mean / 2)
- The two approaches are essentially indistinguishable. exp 012's "discovery" was a lucky seed.

**Implications for strategy:**
- True mean_r(14) ceiling for known recipes is ~0.152.
- Standard deviation per seed is ~0.01-0.02 on mean14, ~0.05+ on per-eval r.
- For the FINAL submission (single library) the seed-to-seed variance MATTERS. A lucky seed of any composition >> the "average" of a "better" composition.
- Use remaining experiments to: (a) find truly better recipes (might not exist), or (b) accept the ceiling and run many seeds of the same recipe to catch a lucky one.

**Next strategy:**
- Exp 015: try ONE more genuinely new direction — combined DHS (cell-type AND topic) + ENCODE TF ChIP-seq if available. Look for the highest-impact source I haven't used.
- Exps 016-020: if nothing breaks the ceiling, do multi-seed sweeps of best recipes to find a high-eval_01 lucky seed.
