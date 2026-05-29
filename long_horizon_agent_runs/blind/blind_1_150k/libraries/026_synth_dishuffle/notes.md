# 026 — Dinucleotide-shuffled cCRE synthetic

## Goal
Test whether composition-matched synthetic (preserves cCRE local
dinucleotide composition but destroys all motifs) substitutes for
motif-embedded random as the diversity component.

## Method
- 67,500 cCRE-fwd + 67,500 different cCRE-RC (135k unique, mixed strand)
- 15,000 dinucleotide-shuffled cCRE windows (Eulerian-walk shuffle preserves
  2-mer composition exactly)
- Total 150k, seed=25

## Result: 0.8895 — within noise of cluster. Eval_08 = 0.9047 (vs 015's 0.9115, **−0.007**).

| eval | 015 (motif synth) | 026 (di-shuf synth) | Δ |
|------|-------------------|---------------------|---|
| 07 | 0.9063 | 0.9050 | −0.001 |
| 08 | 0.9115 | 0.9047 | **−0.007** |
| 13 | 0.9049 | 0.9069 | +0.002 |
| 11 | 0.8187 | 0.8168 | −0.002 |
| 12 | 0.8011 | 0.8003 | −0.001 |

## Key observations
1. **Mean within noise** of motif-embedded synthetic (0.8895 vs 0.8905).
2. **Eval_08 specifically drops 0.007** with di-shuffled synthetic — confirms
   motif content in synthetic IS doing partial work for the eval_08 bonus.
   It's not ONLY composition contrast.
3. **The eval_08 mechanism is dual:**
   - composition contrast (random ACGT or shuffled cCRE both differ from cCRE k-mer
     spectrum) → some bonus
   - motif content (random with embedded JASPAR motifs) → more bonus
4. **Eval_13 slightly UP** with di-shuffled (+0.002). Suggests eval_13 benefits
   from "more cCRE-like" content even in the synthetic slot.

## Theory update (v24 → v25)
The eval_08 super-additivity = (composition contrast) + (motif presence).
About 70 % from composition (broader "this looks different from cCRE")
and 30 % from explicit motifs (~0.007 of 0.024 in 015's case).

This is the most interesting mechanistic finding in the synthetic
exploration.

## Implication
For final library design, use motif-embedded random as synthetic (not
di-shuffled cCRE) — gives full eval_08 bonus. The di-shuffled synthetic
loses 0.007 on eval_08 with negligible compensating gains.

## Next
EXPERIMENT 027 = combine best ingredients for a final design attempt.
