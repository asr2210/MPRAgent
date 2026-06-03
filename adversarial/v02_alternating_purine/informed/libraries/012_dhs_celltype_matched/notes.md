# 012_dhs_celltype_matched — NEW BEST: mean_r=0.162, eval_01=0.174

**Composition (50k):**
- 8k DHS Neural (SK-N-SH match)
- 6k DHS Myeloid/erythroid (K562 match)
- 6k DHS Digestive (HepG2 match)
- 5k cCRE uniform
- 5k Table_S2
- 20k paired flanks

**Result:** mean_r(14) = **0.1620**, eval_01 = **0.1741**. New personal best.

**Vs recipe-baseline (multi-source ~0.151):** +0.011 mean14, +0.015 eval_01.

**Per-eval (selected):**
| eval | exp 7/11 baseline | exp 012 | Δ |
|---|---|---|---|
| eval_01 (chr-GT) | 0.156 | 0.174 | +0.018 |
| eval_06/11 (UKBB K562) | 0.190 | 0.221 | +0.031 |
| eval_03/12 (UKBB) | 0.168 | 0.185 | +0.017 |
| eval_07 (SEI) | 0.16 | 0.162 | +0.00 |
| eval_10 (DHS chr7/13) | 0.16 | 0.123 | -0.037 |
| eval_13 (genomic chr7/13) | 0.13 | 0.129 | 0.00 |

**Cell-type r per eval (mean across UKBB family):**
| | K562 r | HepG2 r | SKNSH r |
|---|---|---|---|
| exp 7/11 baseline | -0.03 to +0.01 | -0.03 to +0.01 | 0.45 |
| exp 012 | +0.03 to +0.08 | +0.03 to +0.08 | 0.45-0.51 |

K562 and HepG2 r LIFTED across all variant-MPRA evals. SKNSH r unchanged or slightly higher. This is exactly what the cell-type-matched DHS hypothesis predicted.

**Theory v10 — cell-type-matched DHS enrichment works:**
The three model output heads (K562/HepG2/SKNSH) are tied to specific cell types. Including DHS sequences pre-labeled as active in matching chromatin programs (Myeloid/erythroid for K562, Digestive for HepG2, Neural for SKNSH) directly trains each head on relevant signal. This is more effective than:
- Asymmetric cCRE stratification (CTCF/DNH3/PLS classes are cell-type-agnostic)
- Table_S2 K562-top-quartile (might be too narrow, training distribution different)

**Trade-off:** eval_10 (DHS chr7/13) dropped 0.04. This is because the topic-weighted DHS that previously dominated eval_10 was replaced with cell-type-specific DHS. Reasonable trade for +0.011 net mean14.

**Next:** verify reproducibility with seed=1 OR push further with bigger DHS-component shares to see how high mean14 can go.
