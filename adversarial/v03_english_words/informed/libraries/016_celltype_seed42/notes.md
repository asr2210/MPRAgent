# 016 — Exp 008 design re-run with SEED=42 (noise floor estimate)

## Plan
Exact copy of exp 008 (best, 0.4283), with only one change: SEED=42 instead
of SEED=0. Measures library randomization noise.

## Result
**eval_01 = 0.4222.** vs exp 008 (SEED=0) = 0.4283. DIFFERENCE: -0.006
from changing only the seed.

K562: 0.593, HepG2: 0.623, SK-N-SH: 0.050. eval_07 = 0.4314 (highest seen
across all experiments, SK-N-SH = 0.065 here breaking 0.06 ceiling).

## What this teaches — MAJOR finding
**Library noise floor ≈ ±0.005-0.01 on mean_r.** Differences smaller than
this between experiments cannot be trusted as signal.

Re-examining recent results through noise lens:
| Exp | mean_r | vs 008 | Above noise? |
|-----|--------|--------|--------------|
| 008 | 0.4283 |  ref   |              |
| 009 | 0.4259 | -0.002 | NO           |
| 010 | 0.4217 | -0.007 | borderline   |
| 011 | 0.4169 | -0.011 | YES (worse)  |
| 012 | 0.4196 | -0.009 | borderline   |
| 013 | 0.4220 | -0.006 | NO           |
| 014 | 0.4184 | -0.010 | borderline   |
| 015 | 0.4220 | -0.006 | NO           |
| 016 | 0.4222 | -0.006 | NOISE ref    |

**Conclusion:** Exps 009, 013, 015 are within noise of exp 008 — meaning
density (5 vs 3), focused pool, and variable density all yield essentially
equivalent libraries to the 008 baseline. Only exp 011 (consensus) is
clearly worse — confirming that PFM stochasticity is critically helpful.

Exp 008 was at the upper end of its noise envelope. True mean_r for the
"3 motifs / 289 cell-type pool / random backbone" family is approximately
0.422-0.428. Single-seed comparisons of similar designs are unreliable.

## Theory T12
- The 008-family plateau is **~0.425 ± 0.005** with this pipeline.
- Six "regression" experiments were actually within noise of equivalent.
- To beat this plateau requires substantively different library structure,
  not parameter tweaks of the same family.

## Implications for next experiments
- Stop fine-tuning 008-family designs.
- Try qualitatively different approaches:
  1. Realistic enhancer SYNTAX (paired motifs with biological spacing)
  2. Much LARGER motif pool (full JASPAR) with composition-aware sampling
  3. Multi-library MIXTURE (motif + DHS + random)
  4. Sequence SCORING + selection: generate 200k candidates, keep best 50k
  5. Promoter-motif structures (TATA + initiator + downstream)

Multi-seed evaluation should be used to validate any candidate that beats
the plateau by less than 0.015.
