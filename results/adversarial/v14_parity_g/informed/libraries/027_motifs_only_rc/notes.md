# 027_ccre_only_rc

Ablation: only cCRE (5 classes × 5k = 25k unique parents) + RC = 50k.
No synthetic motifs, no Sharpr. Pairs with 028 (Sharpr-only+RC) and
019 (motif-only+RC) to isolate per-source contributions to flagship.

Result: eval_01 = -0.0005. Noise.

Triangulation: 019 (synth+RC) = +0.0025, 028 (Sharpr+RC) = +0.0003,
027 (cCRE+RC) = -0.0005. All within noise. No single source dominates.
