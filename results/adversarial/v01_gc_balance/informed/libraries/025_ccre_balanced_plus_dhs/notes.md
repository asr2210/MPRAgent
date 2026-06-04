# 025_ccre_balanced_plus_dhs — notes

**Design**: 25k cCRE class-balanced (~3125 per of 8 classes) + 25k DHS no-label uniform, de-duplicated by genomic midpoint.

**Result**: eval_01 = 0.6858 (vs exp 024 = 0.6921, **-0.006**). Mean = 0.6393 (vs 0.6460, -0.007). Slight decline.

**Comparison**: between exp 022 (0.6827 all-cCRE natural) and exp 024 (0.6921 all-cCRE balanced).

The 25k DHS half pulls the average down — it's a worse source than balanced-cCRE on its own (DHS = 0.6752 < cCRE-balanced = 0.6921). Mixing in 50% of a weaker source predictably regresses partway toward it.

**Confirms H8/H10**: when one source is strictly better, mixing dilutes. Cross-source diversity matters less than within-source structural diversity (class balance).

**Best library remains exp 024 (cCRE class-balanced, 0.6921).**

**Next**: stop trying to add DHS — focus on PUSHING cCRE class-balance further (rare-class oversample, chrom-balance).
