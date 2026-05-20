# 012 — RC-augmented cCRE + motif-embedded (NEW BEST 0.8905)

## Goal
Test whether reverse-complement augmentation gives the model strand-aware
training signal that beats 006's design at the same total library size.

## Method
- 67,500 unique cCRE-centered 200bp windows (forward), seed=11
- 67,500 reverse complements of the same 67,500 cCREs
- 15,000 motif-embedded synthetic (2 JASPAR consensus motifs/seq), seed=11
- Concat, shuffle, write 150k.

## Result: NEW BEST. Mean = 0.8905 (vs 009 0.8888, +0.002; vs 006 0.8883, +0.002).

| eval | 006 (135k unique + 15k motif) | 012 (67.5k×2 + 15k motif) | Δ |
|------|-------------------------------|---------------------------|---|
| 01 | 0.8295 | 0.8308 | +0.001 |
| 02 | 0.9280 | 0.9312 | +0.003 |
| 03 | 0.9207 | 0.9240 | +0.003 |
| 04 | 0.8661 | 0.8686 | +0.003 |
| 05 | 0.8293 | 0.8305 | +0.001 |
| 06 | 0.9284 | 0.9313 | +0.003 |
| 07 | 0.9062 | 0.9061 | 0.000 |
| 08 | 0.9088 | 0.9115 | +0.003 |
| 09 | 0.9452 | 0.9482 | +0.003 |
| 10 | 0.9287 | 0.9311 | +0.002 |
| 11 | 0.8155 | 0.8163 | +0.001 |
| 12 | 0.7984 | **0.8000** | +0.002 |
| 13 | 0.9022 | 0.9052 | +0.003 |
| 14 | 0.9285 | 0.9317 | +0.003 |

## Key observations
1. **RC augmentation is a real lever.** +0.002 mean despite using half the
   unique cCREs (67.5k vs 135k in 006). The 2x effective per-cCRE coverage
   beats the diversity loss.
2. **Wins on 13 of 14 evals.** Consistent improvement across the board.
3. **Eval_12 finally crosses 0.80** for the first time across all
   experiments. Hard evals still hard but trending up.
4. **Eval_08 stays high** (0.9115). The synthetic component still does its
   eval_08 job. The RC of cCREs doesn't dilute the eval_08 bonus.
5. **The prepare.py model is NOT fully strand-invariant.** Augmenting with
   RC gives meaningful new training signal.

## Theory update (v11 → v12)
- Three independent levers now identified:
    1. Function enrichment (random → genome → cCRE): saturated at cCRE,
       gave +0.066 vs random
    2. Synthetic diversity (10% motif-embedded): +0.002 mean, +0.017 eval_08
    3. RC augmentation: +0.002 mean, +0.002 across the board
- They appear roughly **additive**: 003 (0.886) → 006 (+0.002 synthetic) →
  012 (+0.002 RC + synthetic) = 0.890. ~+0.004 total above pure cCRE.
- Hard evals (11, 12) finally moving — eval_12 crossed 0.80.
- The model's strand-aware imperfection is exploitable. Library design
  should account for ML-model properties (architecture, invariances) not
  just biological properties of sequences.

## Next
Two questions to answer:
- **(a) Is synthetic still needed once we have RC?** Test pure 75k×2 RC
  (no synthetic). If ≈ 0.890, synthetic is redundant with RC. If < 0.890,
  the two are independent contributions.
- **(b) Can the 90/10 ratio be re-tuned for RC libraries?** E.g., 60k×2 +
  30k synthetic.

Going with **(a)** first — it tells me whether to focus on RC (only) or
RC + synthetic going forward. Cheap, decisive, and isolates the RC
contribution from the synthetic contribution.

EXPERIMENT 013 = 75,000 unique cCREs × 2 strands = 150,000 (no synthetic).
