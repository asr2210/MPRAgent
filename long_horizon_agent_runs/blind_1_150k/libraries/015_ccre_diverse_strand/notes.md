# 015 — Diverse cCREs with mixed strand assignment

## Goal
Isolate the RC contribution mechanism. Pair-wise (same cCRE in both strands)
vs strand-mixed (different cCREs in each strand) — which is the active
ingredient?

## Method
- 67,500 cCREs (forward only), seed=14
- 67,500 different cCREs (reverse complement), seed=14
- 15,000 motif-embedded synthetic
- Total: 135,000 UNIQUE cCREs (mixed strand) + 15k synthetic = 150k.

## Result: TIED with 012 at 0.8905.

| eval | 006 (135k fwd) | 012 (67.5k × 2) | 015 (135k mixed) | Δ(015−012) |
|------|----------------|------------------|-------------------|------------|
| 01   | 0.8295         | 0.8308           | 0.8330            | +0.002 |
| 02   | 0.9280         | 0.9312           | 0.9303            | −0.001 |
| 03   | 0.9207         | 0.9240           | 0.9230            | −0.001 |
| 04   | 0.8661         | 0.8686           | 0.8657            | −0.003 |
| 05   | 0.8293         | 0.8305           | 0.8328            | +0.002 |
| 06   | 0.9284         | 0.9313           | 0.9307            | −0.001 |
| 07   | 0.9062         | 0.9061           | 0.9063            | 0.000 |
| 08   | 0.9088         | 0.9115           | 0.9115            | 0.000 |
| 09   | 0.9452         | 0.9482           | 0.9467            | −0.002 |
| 10   | 0.9287         | 0.9311           | 0.9313            | 0.000 |
| 11   | 0.8155         | 0.8163           | **0.8187**        | **+0.002** |
| 12   | 0.7984         | 0.8000           | **0.8011**        | +0.001 |
| 13   | 0.9022         | 0.9052           | 0.9049            | 0.000 |
| 14   | 0.9285         | 0.9317           | 0.9308            | −0.001 |

Mean: **0.8905** (012) vs **0.8905** (015) — IDENTICAL to 4 digits.

## Key observations
1. **Pair-wise RC and mixed-strand RC are mean-equivalent.** Both score
   0.8905. The RC mechanism is "model trains on both strands", not "model
   sees paired examples".
2. **Mixed-strand diverse cCREs (015) lifts hard evals**:
   - eval_11: 0.8163 → 0.8187 (+0.002, ALL-TIME HIGH)
   - eval_12: 0.8000 → 0.8011 (+0.001, ALL-TIME HIGH)
3. **Pair-wise (012) has slightly better easy evals** (09, 02): +0.001 to
   +0.002. Maybe pair-wise gives marginally tighter strand learning on
   "easy" sequences.
4. **Net wash on mean.** The two designs trade off slightly across evals
   but average to the same number.

## Theory update (v14 → v15)
- The RC lever is about **strand variation in the training distribution**,
  not specifically about paired examples. Implication: I should use the
  design that maximises *unique cCREs* if it doesn't cost mean — and 015
  has more unique cCREs (135k vs 67.5k), so 015 is the better design
  philosophically (more biological diversity per sequence).
- Hard evals benefit slightly from more unique cCRE diversity. Direction
  confirmed: bigger cCRE pool helps hard evals (incrementally).

## Next
The "cCRE-derived + 15k synthetic" family has plateaued at 0.890. To break
through, I need either:
- A NEW data source (FANTOM5, DNase peaks, conservation)
- A different lever (designed synthetic, k-mer diversity max)
- Different cCRE selection (cell-type, conservation, larger cCREs)

Going with FANTOM5 enhancers — they're a different annotation methodology
(CAGE bidirectional transcription) and could capture regulatory grammar
ENCODE cCREs miss.

EXPERIMENT 016 = 67.5k cCRE-fwd + 67.5k cCRE-RC (different) + 7.5k FANTOM5
enhancers + 7.5k motif = 150k. Tests multi-source biology + synthetic.
