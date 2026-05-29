# 020 — Wide cCREs (≥300bp) only

## Goal
Test if wider cCREs (≥300bp width — presumably more complex regulatory
regions with richer architecture per sequence) outperform the natural
width distribution.

## Method
- Filter to cCREs with end−start ≥ 300bp (488k of 1.06M total cCREs).
- 67,500 wide-cCRE-fwd + 67,500 different wide-cCRE-RC (135k unique).
- 15k motif-embedded synthetic.

## Result: regression. Mean = 0.8880 (vs 015 0.8905, −0.002).

| eval | 015 | 020 | Δ |
|------|------|------|---|
| 01 | 0.8330 | 0.8307 | −0.002 |
| 07 | 0.9063 | 0.9008 | −0.006 |
| 13 | 0.9049 | 0.9011 | −0.004 |
| 11 | 0.8187 | 0.8165 | −0.002 |
| 12 | 0.8011 | 0.7983 | −0.003 |

## Key observations
1. **Yet another stratification regression.** Width filter joins class
   filter (007), GC stratification (018), and dELS-only (017) as failed
   selection strategies.
2. **Eval_07 and eval_13 lose most** — consistent with prior pattern
   that these evals are sensitive to the natural diversity of cCRE
   regions.

## Theory update (v19 → v20)
**Four independent stratification regressions now confirm**: the natural
cCRE distribution is OPTIMAL. Filter / restrict / rebalance — all hurt.
This is one of the most robust patterns in the project.

The implication for library design: **trust the natural diversity of the
ENCODE cCRE catalogue. Don't try to engineer a "better" subset based on
intrinsic features (class, GC, width).**

## Next
EXPERIMENT 021 = test if 95/5 cCRE/synthetic ratio (slightly more cCREs,
less synthetic) helps. Hybrid ratio sweep that hasn't been done with the
RC family. If 95/5 is better, the optimum on the synthetic-fraction
axis is lower than 10%.
