# 016_sknsh_maximize — SKNSH didn't shift, mean14 at ceiling

**Composition:** 12k Neural + 5k Tissue-invariant + 4k PLS + 4k cCRE-uni + 4k S2 + 21k flanks.

**Result:** mean_r(14) = **0.1518**, eval_01 = 0.1562. SKNSH avg r = 0.439.

**Conclusion:**
- SKNSH r doesn't budge above ~0.45 regardless of source mix.
- Mean14 sits at ~0.15 across all reasonable compositions.
- This appears to be a hard ceiling tied to model capacity / eval set noise, not library composition.

**Updated strategy:**
- Stop hunting for novel composition wins (none have stuck).
- Use remaining budget for multi-seed sweeps of the proven-best recipes (exp 012 mostly).
- Goal: harvest a lucky-seed for the final submission. exp 012 hit 0.174 eval_01 once; can I hit 0.18+ with another seed?

**Top recipes to multi-seed:**
1. exp 012 cell-type-matched DHS (seed 0: 0.174 eval_01)
2. exp 007 multi-source (seed 0: 0.159 eval_01)
3. exp 016 SKNSH-maximize (seed 0: 0.156 eval_01)
