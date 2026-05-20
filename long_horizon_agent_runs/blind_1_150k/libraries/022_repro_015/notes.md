# 022 — Noise / replication test (015 with different seed)

## Goal
Quantify training-pipeline noise. Repeat exp 015's recipe with seed=21
(was seed=14). Same composition, different cCRE sample + different
synthetic generation.

## Method
- 67,500 cCRE-fwd + 67,500 different cCRE-RC + 15k motif (= exp 015 recipe)
- Only seed changed: 14 → 21
- Library generation re-samples cCREs and re-generates synthetic

## Result: 0.8875 vs 015 0.8905, Δ = −0.0030.

| eval | 015 | 022 | Δ |
|------|------|------|---|
| 01 | 0.8330 | 0.8294 | −0.004 |
| 02 | 0.9303 | 0.9281 | −0.002 |
| 03 | 0.9230 | 0.9207 | −0.002 |
| 04 | 0.8657 | 0.8645 | −0.001 |
| 05 | 0.8328 | 0.8291 | −0.004 |
| 06 | 0.9307 | 0.9284 | −0.002 |
| 07 | 0.9063 | 0.9028 | −0.004 |
| 08 | 0.9115 | 0.9068 | −0.005 |
| 09 | 0.9467 | 0.9444 | −0.002 |
| 10 | 0.9313 | 0.9296 | −0.002 |
| 11 | 0.8187 | 0.8152 | −0.003 |
| 12 | 0.8011 | 0.7976 | −0.003 |
| 13 | 0.9049 | 0.8997 | −0.005 |
| 14 | 0.9308 | 0.9286 | −0.002 |

## Key observations
1. **Noise floor on mean ≈ ±0.003**. Two runs of the SAME recipe yielded
   means differing by 0.003.
2. **Noise per-eval ≈ ±0.003-0.005.** Individual evals can shift by up
   to half a percentage point just from re-sampling.
3. **This recalibrates ALL prior 'wins'**:
   - 003 (cCRE) → 006 (90/10 hybrid): +0.002 — WITHIN NOISE, not a confirmed win
   - 006 → 012 (RC): +0.002 — WITHIN NOISE, not a confirmed win
   - 012 → 015 (mixed strand): 0.000 — confirmed equivalent
   - 015 → 016 (FANTOM5 mix): +0.0002 — WITHIN NOISE
   - 009 (dense motif): +0.0005 — WITHIN NOISE
4. **What IS robustly above noise (Δ ≥ 0.005)**:
   - random (0.820) → genome (0.873): +0.053 ✓ HUGE
   - genome → cCRE (0.886): +0.013 ✓ real
   - cCRE → cCRE+10%synth: +0.002 (within noise but consistent)
   - All 4 cCRE-stratification regressions: −0.004 to −0.015 ✓ real
   - 50/50 hybrid eval_08 super-additivity: +0.024 on eval_08 ✓ real

## Theory update (v21 → v22)
**The 0.886-0.892 'ceiling' is fuzzy.** What I thought were distinguishable
designs (012, 015, 016, 019 all reported 0.889-0.891) are all
indistinguishable within training noise. The true performance of the
cCRE + RC + 10% synthetic family is **0.889 ± 0.003**.

This is a CRUCIAL meta-finding: many marginal "improvements" in the second
half of the project were noise. The genuine library design wisdom is:
- Use cCRE-derived bulk (huge over random/genome)
- Don't filter cCRE distribution (each stratification regression confirmed)
- Add ~10% synthetic (modest but consistent direction)
- RC augmentation is probably neutral-to-small-positive (within noise)
- Synthetic source is fungible

## Next
Given noise is ±0.003, I should:
- STOP fine-tuning the cCRE + RC + 10% synthetic recipe; it's saturated
- Try genuinely new levers (cell-type, multi-source bulk, conservation)
- Take a few replicates of the best design to nail down its mean
- Build a final "best known" library

EXPERIMENT 023 = KITCHEN-SINK combining all winning ingredients to see if
multi-lever stacking gives ANY compound benefit above noise. If not,
remaining experiments will go to replicates and final library design.
