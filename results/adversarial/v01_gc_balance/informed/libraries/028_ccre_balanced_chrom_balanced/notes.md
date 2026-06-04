# 028_ccre_balanced_chrom_balanced — notes

**Design**: cCRE 8-class balanced (exp 024 design) PLUS per-chromosome balanced within each class. 8 classes × 23 chroms (autosomes + chrX) × ~272 sequences/cell = 50,000. Drops chrY/chrM (very sparse cCRE coverage).

**Result**: eval_01 = **0.6940** (vs exp 024 = 0.6921, **+0.002**). Mean = 0.6479 (vs 0.6460, +0.002).

🎯 **NEW BEST** — chrom-balance is orthogonal to class-balance and lifts.

Per-eval vs exp 024:
| eval | exp 024 | exp 028 | Δ |
|---|---|---|---|
| eval_01 | 0.6921 | 0.6940 | +0.002 |
| eval_04/09 | 0.5966 | 0.6072 | **+0.011** |
| eval_07 | 0.7576 | 0.7548 | -0.003 |
| eval_13 | 0.7479 | 0.7444 | -0.004 |

The biggest lift is again on eval_04/09 (hardest evals). Chrom-balance pulls in more chrX and small-chromosome regions, which were under-sampled in natural cCRE density. Slight regression on eval_07/13.

**H13**: Genomic distribution of training sequences matters for generalization. Uniform sampling within class still over-represents large chromosomes; explicitly balancing across chromosomes (within each class) brings in regulatory grammar from gene-poor / late-replicating regions that the class-balanced library was missing.

**Best library: exp 028 (cCRE 8-class × 23-chrom balanced, 0.6940).**

**Next**:
- Exp 029: 7-class chrom-balanced (no dELS, combining the wins from 026 and 028) — does dropping dELS still help on top of chrom-balance?
- Exp 030: multi-seed verification of best recipe.
