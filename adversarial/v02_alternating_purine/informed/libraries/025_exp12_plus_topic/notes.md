# 025_exp12_plus_topic — topic addition fails to recover eval_10

**Composition:** 7k Neural + 5k Myeloid + 5k Digestive + **3k DHS-topic (NMF-weighted)** + 5k cCRE-uni + 5k S2 + 20k flanks.

**Result:** mean14 = **0.1489**, eval_01 = 0.1538, eval_10 = 0.1239 (NOT recovered).

**Hypothesis falsified:** Adding 3k NMF-topic-weighted DHS did not recover eval_10. The exp 012 family's eval_10 drop is structural (eval_10 = chr7/13 DHS), not source-mix-driven. Topic-weighted DHS likely already overlaps with Neural/Myeloid/Digestive components, so we added redundancy, not coverage.

**eval_10 numbers across recipes:**
- exp 001 DHS-topic-only: 0.169
- exp 002 DHS-topic + flanks: 0.186
- exp 012 cell-type DHS: 0.123
- exp 025 (this) cell-type + topic: 0.124

So 3k topic was insufficient; topic-dominated recipes (exp 001/002) hit 0.17-0.19 on eval_10 but lost on others. Trading 3k of cell-type DHS for topic gave no eval_10 lift because there isn't enough topic to dominate.

**Conclusion:** Stop trying to fix eval_10 inside cell-type recipe. The eval_10 penalty (~0.04) is a structural cost of cell-type-matched DHS. Net mean14 is still better than topic-pure recipes.

**Next:** 4-5 experiments left. Keep sweeping exp 012 seeds for lucky-seed harvest. Current best: 019_exp12_seed4 (mean14=0.164, eval_01=0.177).
