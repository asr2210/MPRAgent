# 010_s2_k562_actives — K562-active anchors via Table_S2 stratification

**Composition:** 10k uniform + 4k CTCF + 4k DNH3 + 3k S2 K562-top-quartile + 2k S2 random + 4k DHS + 3k synth + 20k flanks (mirrors exp 007 with S2 split into K562-top vs random)

**Result:** mean_r(14) = **0.1361**, eval_01 = 0.1361. *Regression vs exp 007 (0.1532, 0.1592).*

**The disturbing finding:** This recipe is STRUCTURALLY ALMOST IDENTICAL to exp 007. Only difference: 5k random Table_S2 → 3k K562-top-quartile + 2k random S2 (and a different sampling order due to seen-set state). Yet:
- K562 r flipped from +0.07 to -0.02 on eval_06/11
- eval_06/11 mean_r dropped from 0.228 to 0.156

This means EITHER:
(a) K562-top-quartile selection is actively harmful (selects atypical sequences), OR
(b) Exp 007's eval_06/11 = 0.228 was a high-variance lucky-seed result, not a reproducible recipe-level effect

**Looking back at K562 r history (eval_06):**
- Exp 1: -0.032 | Exp 2: -0.039 | Exp 3: -0.028 | Exp 4: -0.018 | Exp 5: -0.053
- Exp 6: -0.018 | **Exp 7: +0.071 (OUTLIER)** | Exp 8: +0.026 | Exp 9: -0.039 | Exp 10: -0.021

Exp 7's K562=+0.071 is 4σ above the typical ~-0.03. The harness probably has substantial seed-to-seed variance, and that lucky shot inflated our "best" estimate.

**This is the most important finding so far.** I've been chasing a high-variance signal. My subsequent comparisons are unreliable to within ±0.04 on individual eval cell-type r.

**Next experiment:** REPRODUCIBILITY TEST. Re-run exp 007's exact recipe with seed=1. If mean_r ~ 0.14 → exp 007 was lucky; throw out the +0.022 lift. If mean_r ~ 0.15 → exp 007 was real, exp 010 perturbation was the culprit.
