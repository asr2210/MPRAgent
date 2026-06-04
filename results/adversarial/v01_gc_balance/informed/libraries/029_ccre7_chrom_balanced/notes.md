# 029_ccre7_chrom_balanced — notes

**Design**: cCRE 7-class (no dELS) × 23-chrom balanced. Combines exp 026 (drop dELS) and exp 028 (chrom-balance) recipes.

**Result**: eval_01 = 0.6933 (vs exp 028 = 0.6940, **-0.001**). Mean = 0.6477 (vs 0.6479, -0.0002). Essentially flat.

eval_04/09 = 0.6090 (vs 0.6072, +0.002) — marginal lift on hardest.

**Conclusion**: dropping dELS provides no additional lift when chrom-balance is already applied. The dELS contribution at 8.7% (29) or 12.5% (28) or 0% — all flatten out. The "specialized class" hypothesis (H12) is more about RAISING rare classes than about LOWERING dELS.

**Best library remains exp 028 (cCRE 8-class × 23-chrom balanced, 0.6940).**

**Next (exp 030)**: final experiment to verify and push.
