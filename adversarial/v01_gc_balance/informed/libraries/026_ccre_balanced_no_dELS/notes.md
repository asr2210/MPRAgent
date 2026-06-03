# 026_ccre_balanced_no_dELS — notes

**Design**: 7 cCRE classes balanced (~7143 each), dELS dropped entirely. Tests whether dELS at 12.5% (as in exp 024) carried distinct grammar or was redundant with the other 7 classes.

**Result**: eval_01 = **0.6925** (vs exp 024 = 0.6921, **+0.0004**). Mean = 0.6464 (vs 0.6460, +0.0004). **Statistical flat** — within noise.

**Interpretation**: dELS at 12.5% in exp 024 contributed essentially nothing on top of the other 7 classes. pELS is structurally similar to dELS (both enhancer-like) and likely covers most of the same grammar at 12.5% each.

**H12**: Beyond the broad coverage of distinct classes, additional fine-class diversity within enhancer-like grammars (dELS vs pELS) is redundant. The 5 "specialized" classes (PLS, CA-TF, CA-H3K4me3, TF, CA-CTCF) are the ones that lift; the 3 enhancer-like classes (dELS, pELS, CA which is "accessible-only") are mutually redundant.

🎯 **NEW BEST (tied): exp 026 (0.6925) and exp 024 (0.6921)** — both are within noise. Exp 026 is structurally simpler (drops the largest natural class without loss).

**Next**: test other ablations.
- Exp 027: drop pELS too (6 classes, ~8333 each) — does removing both enhancer-like classes hurt?
- Exp 028: cCRE balanced + per-class oversample (e.g., 3x for rarest classes via with-replacement) — does going further than "balanced" toward rare classes help?
- Exp 029: chrom-balanced within class — combine class-balance with chromosome-balance
