# 013_dhs_celltype_aggressive — scale-up regression

**Composition:** 10k Neural + 8k Myeloid + 8k Digestive + 4k cCRE-uni + 0k S2 + 20k flanks.

**Result:** mean_r(14) = **0.1512**, eval_01 = 0.1539. *Regression vs exp 012's 0.162.*

**Key deltas vs exp 012:**
- eval_06/11 K562 r: 0.078 → 0.019 (lost the K562 lift)
- eval_06/11 mean: 0.221 → 0.181
- eval_07: 0.162 → 0.174 (interesting — gained without Table_S2)
- eval_13: 0.129 → 0.161 (better, perhaps fewer cCRE = more genomic flank diversity)

**Interpretation:**
Dropping 5k Table_S2 cost ~0.04 K562 r on UKBB evals. Adding more cell-type DHS (6k extra) didn't compensate. **Theory v11 update:** Table_S2 is itself a K562 trainer (it has measured K562 activity values).

Wait — but exp 005 (50k pure Table_S2) had eval_06 = 0.106 and eval_11 = 0.106. So Table_S2 ALONE doesn't give K562 r. The 0.078 K562 r in exp 012 came from the COMBINATION of cell-type DHS + Table_S2 + cCRE uniform. Dropping any of them disrupts the balance.

**Possible explanations:**
1. **Interaction-based:** all three sources contribute orthogonal K562-relevant features; need all of them.
2. **Lucky seed:** exp 012's 0.078 K562 r might be at the upper end of a noisy distribution.
3. **Diversity ceiling:** too much DHS (overfit to chromatin-program signature) collapses representation; need source diversity.

**Next:**
- Exp 014: Verify exp 012 reproducibility with seed=1.
- If reproducible: tweak exp 012 by ADDING (not replacing) — e.g., extend cell-type DHS while keeping S2 by reducing flanks.
- If not reproducible: revert to the "recipe baseline ~0.151" view and try a totally different angle.
