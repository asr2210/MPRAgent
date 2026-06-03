# 013_dhs_gc_extreme — notes

**Design**: pre-sample 200k DHSs from non-label-aligned pool, compute GC content per 200bp window, take 25k lowest + 25k highest GC DHSs. Total 50k from extreme tails (GC <0.36 or >0.57; pool natural GC = 0.46±0.10).

**Result**: eval_01 = 0.6384 (vs exp 005 = 0.6752, **-0.037**). Mean = 0.5882 (vs 0.6282, -0.040).

**Counter-test for H5 — confirmed**. Concentration on a totally different axis (GC composition) hurts by similar magnitude as motif-density concentration (exp 011). Generalizing across 4 concentration experiments now (007 numsamples, 008 forced balance, 011 motif-dense, 013 GC-extreme):

| concentration axis | eval_01 | Δ vs exp 005 |
|---|---|---|
| numsamples-weighted | 0.6185 | -0.057 |
| topic-balanced | 0.6631 | -0.012 |
| motif-dense | 0.6084 | -0.067 |
| GC-extreme | 0.6384 | -0.037 |
| **uniform (exp 005)** | **0.6752** | — |

**H5 confidently confirmed**: any selection axis that concentrates the library reduces compositional/contextual variance and hurts. The model relies on the natural distribution of regulatory sequences.

**Next direction**: stay with uniform sampling, vary the FILTER (which DHSs to exclude) rather than the SELECT (which to favor). Filter is set-membership; select is rank/weight.
