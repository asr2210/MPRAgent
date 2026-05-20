# 014 — Partial RC hybrid (90k unique cCREs, half RC-paired)

## Goal
Test the middle of the unique-vs-RC trade-off:
- 006 (135k unique fwd + 15k motif): 0.8883 — all unique, no RC
- 012 (67.5k × 2 + 15k motif):        0.8905 — half unique, full RC
- 014 (45k × 2 + 45k fwd + 15k motif): ?    — 90k unique, partial RC

## Result: 0.8879 — slightly below 006, well below 012.

| eval | 012 | 014 | Δ(014−012) |
|------|------|------|------------|
| 01 | 0.8308 | 0.8307 | 0.000 |
| 02 | 0.9312 | 0.9282 | −0.003 |
| 03 | 0.9240 | 0.9206 | −0.003 |
| 04 | 0.8686 | 0.8639 | −0.005 |
| 05 | 0.8305 | 0.8304 | 0.000 |
| 06 | 0.9313 | 0.9285 | −0.003 |
| 07 | 0.9061 | 0.9024 | −0.004 |
| 08 | 0.9115 | 0.9085 | −0.003 |
| 09 | 0.9482 | 0.9440 | −0.004 |
| 10 | 0.9311 | 0.9283 | −0.003 |
| 11 | 0.8163 | 0.8164 | 0.000 |
| 12 | 0.8000 | 0.7983 | −0.002 |
| 13 | 0.9052 | 0.9021 | −0.003 |
| 14 | 0.9317 | 0.9288 | −0.003 |

## Key observations
1. **Full RC pairing (012) > partial RC (014) on 11 evals.** Tied on
   3. Never won. The model benefits more from pairing every cCRE on both
   strands than from a mix of "some paired, some not".
2. **014 mean ≈ 006 mean** (0.8879 vs 0.8883). Adding partial RC at the
   cost of cCRE diversity is approximately a wash.
3. **Hard evals 11/12 plateau** (0.816, 0.798).
4. The trade-off curve has a clear local optimum at 012's design point
   (all unique cCREs RC-paired).

## Theory update (v13 → v14)
- Within the "cCRE + RC + synthetic" family, exp 012's full-pairing design
  is the leader. Partial pairing dilutes the value of both:
    - More unique cCREs but fewer RC pairings → loses some of the RC lift
    - Fewer unique cCREs in pairs but more total → loses cCRE coverage
  Neither half-measure recovers the full lift.
- The RC effect likely depends on the model seeing BOTH orientations of
  the SAME sequence — it can compare and learn strand-symmetric grammar
  more effectively than from disparate strand-mixed examples. Confirmed
  by partial-pairing regression.

## Next
Going to test: does **cCRE diversity beat strand-pairing** if we use 135k
unique cCREs but with each one in a random orientation (forward or RC)?
Same total cCRE-derived count as 012 (135k), more unique cCREs (135k vs
67.5k), but no pair-wise strand training.

EXPERIMENT 015 = 67.5k cCREs (forward) + 67.5k DIFFERENT cCREs (RC) +
15k motif = 150k. 135k unique, mixed strand, with synthetic.

If 015 > 012: cCRE diversity matters more than pair-wise strand training
If 015 < 012: pair-wise strand training is the active ingredient
