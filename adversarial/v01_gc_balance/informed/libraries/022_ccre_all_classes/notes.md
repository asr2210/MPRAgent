# 022_ccre_all_classes — notes

**Design**: 50k uniform sample from ALL 8 cCRE classes (2.35M total): dELS 62.5%, pELS 10.6%, CA 10.5%, CA-CTCF 5.4%, TF 4.5%, CA-H3K4me3 3.4%, PLS 2.0%, CA-TF 1.1%. Window centered on midpoint, ±100bp.

**Result**: eval_01 = **0.6827** (vs exp 005 = 0.6752, **+0.0075**). Mean = 0.6381 (vs 0.6282, +0.010).

🎯 **First experiment to beat exp 005.** Both eval_01 and mean improve.

Per-eval lift over exp 005:
| eval | exp 005 | exp 022 | Δ |
|---|---|---|---|
| eval_01 | 0.6752 | 0.6827 | +0.008 |
| eval_04/09 | 0.5374 | 0.5704 | +0.033 |
| eval_07 | 0.7597 | 0.7554 | -0.004 |
| eval_13 | 0.7509 | 0.7451 | -0.006 |

The lift comes from broad gains (eval_01/02/03/04/05/06/10/11/12/14), with a small cost on eval_07/13.

**Comparison to exp 010 (cCRE dELS only) = 0.6671**:
- All-classes (+0.016 vs dELS-only) — class diversity beyond enhancer-like helps.
- Adding pELS+PLS+CTCF+TF+CA classes contributes regulatory grammar diversity that DHS no-label-topic uniform can't fully capture.

**Theory update H9 → H10**:
> *Diversity of regulatory program TYPES* (enhancer + promoter + CTCF + TF + accessible-only) beats either pure-enhancer or DHS topic-exclusion. cCRE integrates multiple chromatin marks (DNase + H3K27ac + H3K4me3 + CTCF + TF) and the integration is more informative than DNase-only DHSs filtered by NMF.

**Next**: push this further.
- Exp 023: cCRE class-balanced (forced equal per class) — tests whether the natural class proportions are optimal or if explicit balancing helps.
- Exp 024: cCRE all-classes coordinate-intersected with no-label DHS — combine both filters.
- Exp 025: cCRE all-classes EXCLUDING PLS (drop core promoters which may be label-aligned).
