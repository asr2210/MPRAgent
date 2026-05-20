# 009 — 90/10 cCRE + dense motif-embedded (5 motifs/seq)

## Goal
Refine the winning hybrid (exp 006) by increasing the motif density in the
synthetic component from 2 → 5 motifs per 200bp.

## Method
- 135,000 cCRE-centered windows (same pool & filter), seed=8
- 15,000 dense motif-embedded synthetic: random ACGT background + 5 JASPAR
  consensus motifs at non-overlapping positions (with retry up to 50; if
  no free slot found, stamps over), seed=8
- Concatenated, shuffled, written.

## Result: marginal new best. Mean = 0.8888 (vs 006 0.8883, +0.0005).
Runtime 2718s.

| eval | 006 (2 motifs) | 009 (5 motifs) | Δ |
|------|----------------|----------------|---|
| 01 | 0.8295 | 0.8310 | +0.002 |
| 02 | 0.9280 | 0.9293 | +0.001 |
| 03 | 0.9207 | 0.9223 | +0.002 |
| 04 | 0.8661 | 0.8647 | −0.001 |
| 05 | 0.8293 | 0.8308 | +0.002 |
| 06 | 0.9284 | 0.9296 | +0.001 |
| 07 | 0.9062 | 0.9058 | 0.000 |
| 08 | 0.9088 | 0.9062 | −0.003 |
| 09 | 0.9452 | 0.9453 | 0.000 |
| 10 | 0.9287 | 0.9295 | +0.001 |
| 11 | 0.8155 | 0.8166 | +0.001 |
| 12 | 0.7984 | 0.7995 | +0.001 |
| 13 | 0.9022 | 0.9028 | +0.001 |
| 14 | 0.9285 | 0.9298 | +0.001 |

## Key observations
1. **5 motifs barely improves over 2 motifs** (+0.0005 on mean). Most evals
   gain +0.001–+0.002, which is within noise.
2. **Eval_08 actually drops slightly** (−0.003) with more motifs. Consistent
   with the "eval_08 needs random-like backbone" hypothesis — dense
   motifs reduce the random-like character.
3. **Motif density is a flat curve above 2**. Diminishing returns on motif
   density. Stop tuning this knob.
4. **All other evals tiny gains**. The richer motif content gives a small
   training signal boost, but the effect is small.

## Theory update (v8 → v9)
- Motif density in the synthetic component above 2 per 200bp is essentially
  flat. The active ingredient of motif-embedded synthetic is "presence of
  motifs in a random-like backbone", not "as many motifs as possible".
- Eval_08 has a soft trade-off: more motifs ↑ general mean slightly, but
  ↓ eval_08 slightly. Suggests eval_08 specifically benefits from
  more random / less structured backbones.
- **Ceiling of 0.888 holds.** The 90/10 cCRE + synthetic family is
  saturating. Need a fundamentally different lever to break through.

## Next
Going to test a new lever: **cCRE neighbourhoods** (50k unique cCREs ×
3 windows each = 150k). Each cCRE contributes its center window plus
two flanking windows (±200bp shifts). Tests whether sequences *adjacent
to* regulatory elements carry additional generalisable signal beyond
the cCRE centres themselves. Could lift hard evals if they test
extended regulatory context, or could regress if flanks dilute the
high-information centres.
