# 013 — Pure RC-augmented cCRE (no synthetic)

## Goal
Isolate the RC contribution from the synthetic contribution. If pure RC
matches 012 (0.890), synthetic is redundant. If pure RC < 012, synthetic
adds independent value.

## Method
- 75,000 unique cCRE-centered 200bp windows (forward), seed=12
- 75,000 reverse complements of the same 75,000 cCREs
- Total 150k. No synthetic.

## Result: mean = 0.8832. Loses by −0.007 vs 012, and even by −0.003 vs pure cCRE (003).

| eval | 003 cCRE | 012 RC+motif | 013 RC only | Δ(013−012) |
|------|----------|--------------|-------------|------------|
| 01   | 0.8288 | 0.8308 | 0.8265 | −0.004 |
| 02   | 0.9270 | 0.9312 | 0.9245 | −0.007 |
| 03   | 0.9193 | 0.9240 | 0.9162 | −0.008 |
| 04   | 0.8657 | 0.8686 | 0.8629 | −0.006 |
| 05   | 0.8285 | 0.8305 | 0.8264 | −0.004 |
| 06   | 0.9274 | 0.9313 | 0.9249 | −0.006 |
| 07   | 0.9025 | 0.9061 | 0.9000 | −0.006 |
| 08   | 0.8922 | 0.9115 | 0.8876 | **−0.024** |
| 09   | 0.9465 | 0.9482 | 0.9432 | −0.005 |
| 10   | 0.9272 | 0.9311 | 0.9228 | −0.008 |
| 11   | 0.8145 | 0.8163 | 0.8126 | −0.004 |
| 12   | 0.7959 | 0.8000 | 0.7937 | −0.006 |
| 13   | 0.9038 | 0.9052 | 0.8982 | −0.007 |
| 14   | 0.9276 | 0.9317 | 0.9251 | −0.007 |

## Key observations
1. **RC alone underperforms pure cCRE** (−0.003 mean). Showing each of 75k
   cCREs on both strands does NOT compensate for losing half the unique cCRE
   coverage (150k → 75k).
2. **Eval_08 takes the biggest hit** (−0.024). Without synthetic, the
   eval_08 bonus disappears. Confirms the eval_08 effect specifically needs
   synthetic / non-cCRE-distribution sequences.
3. **Pure cCRE > pure RC**: the natural cCRE catalogue with 150k unique
   examples is more informative than 75k unique with both strands. Unique
   cCRE diversity is more valuable per sequence than strand augmentation.
4. **RC is only beneficial WHEN combined with synthetic.** 012 (RC+synth)
   wins; 013 (RC only) loses. This means RC and synthetic are
   non-redundant — they capture different aspects.

## Theory update (v12 → v13)
- The 012 lift (0.890) is the result of TWO mostly-independent
  contributions:
    - Synthetic component: gives eval_08 bonus (+0.017 → +0.024 on that eval)
    - RC augmentation: gives small lift on most other evals (+0.001 to +0.003)
- Neither alone exceeds pure cCRE on mean. Only together do they win.
- This is counter-intuitive: RC alone loses, synthetic alone gives +0.002,
  but RC + synthetic gives +0.004 — super-additive at the mean level.
- Best guess at mechanism: the synthetic component prevents the model from
  over-fitting to cCRE-strand-orientation peculiarities, and RC then
  gives the model strand-symmetric training data to learn strand-aware
  but symmetric grammar. Both ingredients together unlock the full
  generalisation.

## Next
Going to test whether more unique cCREs (less aggressive RC) can push
higher. The trade-off curve has two known points:
- 006: 135k unique forward + 15k motif → 0.888 (lots of unique cCREs, no RC)
- 012: 67.5k unique × 2 + 15k motif → 0.890 (fewer unique, full RC)
- 013: 75k unique × 2 → 0.883 (no synthetic)

Experiment 014 will test the middle: 45k cCREs RC-paired (= 90k seqs) +
45k unique forward-only cCREs + 15k motif = 150k. Gives 90k UNIQUE cCREs
(more than 012) with HALF of them RC-paired. Tests if a middle point
on the unique-vs-RC trade-off wins.
