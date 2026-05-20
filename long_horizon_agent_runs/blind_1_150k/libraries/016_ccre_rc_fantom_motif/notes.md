# 016 — Multi-source: cCRE (mixed strand) + FANTOM5 + motif

## Goal
Test whether FANTOM5 CAGE-defined enhancers (alternative annotation
methodology — bidirectional transcription) add information beyond ENCODE
cCREs when included as part of the diversity component.

## Method
- 67,500 cCREs (forward only), seed=15
- 67,500 different cCREs in reverse complement
- 7,500 FANTOM5 enhancers, 200bp centered, filtered
- 7,500 motif-embedded synthetic
- Total: 135k unique cCREs (mixed strand) + 7.5k FANTOM5 + 7.5k motif = 150k

## Result: marginal new best. Mean = 0.8907 (vs 015 0.8905, +0.0002).

| eval | 015 | 016 | Δ |
|------|------|------|---|
| 01 | 0.8330 | 0.8320 | -0.001 |
| 02 | 0.9303 | 0.9313 | +0.001 |
| 03 | 0.9230 | 0.9244 | +0.001 |
| 04 | 0.8657 | 0.8658 | 0.000 |
| 05 | 0.8328 | 0.8317 | -0.001 |
| 06 | 0.9307 | 0.9316 | +0.001 |
| 07 | 0.9063 | 0.9058 | -0.001 |
| 08 | 0.9115 | 0.9101 | -0.001 |
| 09 | 0.9467 | 0.9469 | 0.000 |
| 10 | 0.9313 | 0.9314 | 0.000 |
| 11 | 0.8187 | 0.8176 | -0.001 |
| 12 | 0.8011 | 0.8012 | 0.000 |
| 13 | 0.9049 | 0.9078 | +0.003 |
| 14 | 0.9308 | 0.9319 | +0.001 |

## Key observations
1. **FANTOM5 ≈ motif-embedded as a diversity component.** Splitting the 15k
   diversity slot between FANTOM5 and motif gives essentially the same
   mean (+0.0002, within noise) as 100% motif.
2. **Eval_13 had the biggest gain** (+0.003). FANTOM5's CAGE-defined
   regulatory regions might particularly suit whatever eval_13 tests.
3. **Hard evals essentially unchanged.** Eval_11 down a hair (-0.001),
   eval_12 flat.
4. **The ceiling at 0.890 is robust.** This is the third design within
   ±0.0002 of 0.890 (012, 015, 016). Suggests genuine plateau.

## Theory update (v15 → v16)
- Diversity component sources (motif-embedded, uniform random, FANTOM5)
  are largely interchangeable. ≈ +0.002 mean over pure cCRE, ≈ +0.020
  on eval_08.
- The 0.890 mean appears to be the ceiling for the
  "cCRE-derived bulk + 15k diversity component" family of designs across
  many variations.
- To break 0.890, need either:
  (a) higher information density per cCRE (better cCRE selection)
  (b) a genuinely new training signal (cell-type-specific, conservation,
      activity-extreme, designed synthetic regulatory grammar)
  (c) algorithmic diversity max within cCRE pool

## Next
Going to test (a): does focusing on dELS-only (the dominant 74% cCRE
class) — which exp 007 implied was the high-information class — give
better signal per sequence?

EXPERIMENT 017 = 67.5k dELS-only cCREs (fwd) + 67.5k different dELS-only
cCREs (RC) + 15k motif. Same recipe as 015 but restricted to dELS class.

If higher: dELS is the active ingredient, focus the library on it.
If equal: all cCRE classes contribute equally per sequence.
If lower: other cCRE classes add value beyond pure dELS.
