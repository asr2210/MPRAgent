# 007 — Class-balanced cCREs

## Goal
Test whether forcing equal representation of PLS / pELS / dELS / CTCF-other
cCRE classes (overriding the natural 74%-dELS distribution) exposes the
model to under-represented promoter and insulator grammars, lifting the
hard evals (11, 12).

## Method
- Bucket cCREs into 4 classes by `class` column.
- Apply same window/filter as exp 003.
- Sample 37,500 per bucket (with replacement when a bucket falls just
  short — PLS-types had 37,331 survivors, Other had 37,216).
- Concat, shuffle, write.

Bucket sizes (after filter):
- PLS-types:  37,331  (37,500 sampled, ~0.5% duplicates)
- pELS-types: 134,608
- dELS-types: 536,243
- Other:      37,216  (37,500 sampled, ~0.8% duplicates)

## Result: REGRESSION. Mean = 0.8820 (vs 0.8862 pure cCRE, **−0.004**).

| eval | natural cCRE | class-balanced | Δ |
|------|--------------|----------------|---|
| 01 | 0.8288 | 0.8251 | −0.004 |
| 02 | 0.9270 | 0.9234 | −0.004 |
| 03 | 0.9193 | 0.9154 | −0.004 |
| 04 | 0.8657 | 0.8657 | 0.000 |
| 05 | 0.8285 | 0.8247 | −0.004 |
| 06 | 0.9274 | 0.9236 | −0.004 |
| 07 | 0.9025 | 0.8918 | **−0.011** |
| 08 | 0.8922 | 0.8946 | +0.002 |
| 09 | 0.9465 | 0.9442 | −0.002 |
| 10 | 0.9272 | 0.9201 | −0.007 |
| 11 | 0.8145 | 0.8107 | −0.004 |
| 12 | 0.7959 | 0.7931 | −0.003 |
| 13 | 0.9038 | 0.8920 | **−0.012** |
| 14 | 0.9276 | 0.9240 | −0.004 |

## Key observations
1. **Losses on 13 of 14 evals**, mostly −0.003 to −0.012. dELS-types are
   evidently the most informative class per sequence; over-representing
   PLS/CTCF (combined 50% of library vs 4% naturally) hurts.
2. **Hard evals stayed hard**: 11 lost a bit (−0.004), 12 lost a bit
   (−0.003). Class balance is NOT the lever for these.
3. **Eval_08 gained slightly** (+0.002). Mild diversity effect — having
   underrepresented classes contributes a small "diversity premium" on
   eval_08, consistent with the broader pattern.
4. **Eval_07 and eval_13 lost most** (−0.011, −0.012). These were exactly
   the evals that benefited most from random→cCRE in exp 003. They
   evidently need many *dELS* (distal enhancer) examples specifically.

## Theory update (v6 → v7)
- The natural cCRE class distribution is GOOD. dELS dominance (74%) is
  serving the model's generalisation, not biasing it. Forcing balance is a
  net loss.
- This is itself informative: it implies that within "regulatory elements",
  distal enhancers are the broadest training signal for the model.
- Hard evals (11, 12) are not unlocked by class diversity. Still need a
  new lever.
- Eval_08 continues to be the diversity-bonus eval — any heterogeneity
  added (random, motif-embedded, class diversity) lifts it slightly.

## Next
Class diversity is dead. The 90/10 hybrid (exp 006) remains the best
design. Next: explore alternative diversity sources or a completely
different lever for the hard evals.

Going with experiment 008 = **90/10 cCRE + non-cCRE genomic tiles**
(replacing motif-embedded as the diversity source with random
non-cCRE 200bp windows from hg38). Direct comparison to exp 006:
biological diversity (genome) vs synthetic diversity (motif-embedded).
If genome > synthetic, biological diversity is preferred. If equal,
diversity source doesn't matter — pick whichever is easier. If
genome < synthetic, the synthetic motif content was specifically
contributing.
