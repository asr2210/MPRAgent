# 027_ccre_specialized_only — notes

**Design**: 5 "specialized" cCRE classes only, balanced at 10k each: PLS (promoters), CA-CTCF (CTCF anchors), TF (TF anchors), CA-H3K4me3 (H3K4me3-marked accessible), CA-TF (TF+accessible). Drops dELS, pELS, CA (the 3 broad enhancer/accessible classes).

**Result**: eval_01 = 0.6917 (vs exp 026 = 0.6925, **-0.001**). Mean = 0.6450 (vs 0.6464, -0.001). **Essentially flat.**

**Surprising**: dropping ALL enhancer-like classes (which make up 75% of cCREs by number) costs almost nothing. The 5 specialized classes alone capture nearly all the predictive grammar.

**H12 confirmed (strong)**: Beyond class-balance, the marginal contribution of dELS/pELS/CA is ~0. The "specialized" classes carry the lift; enhancer-like classes are mutually redundant AND redundant with the specialized ones.

**Implication**: the cCRE class-balance lift over natural proportions is NOT about having all 8 classes — it's about getting the 5 specialized classes to high enough fraction. Natural proportions put them at only ~13% combined; balanced puts them at 62.5%; specialized-only puts them at 100% (with no additional lift).

**Best tied: exp 026 (0.6925), exp 024 (0.6921), exp 027 (0.6917) — within noise.**

**Next**: orthogonal axis. Chrom-balance within class (exp 028) — does forcing uniformity across chromosomes within the class-balanced library add anything? Then multi-seed for the best recipe.
