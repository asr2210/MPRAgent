# 002 — DHS + Dirichlet composition (50/50)

## What I built
50,000 sequences = 25,000 stratified DHS (Meuleman 160k subset) + 25,000 Dirichlet-composition synthetic. Dirichlet sequences: per-sequence base frequencies ~ Dirichlet(1,1,1,1), then 200 bp i.i.d. from that composition.

## Predicted
eval_01 ~ 0.075-0.085 if additive; ~0.075 if not.

## Result
- eval_01 = 0.0765 (exp 001: 0.0739, dirichlet baseline: 0.0768, random: 0.0711) — almost exactly = dirichlet alone, slight win over exp 001.
- eval_08 = 0.0728 (exp 001: 0.0627, dirichlet: 0.0745, random: 0.0605) — biggest jump vs exp 001 (+0.010).
- eval_07 = 0.1464, eval_10 = 0.1300, eval_13 = 0.1443 — biology-sensitive evals still near upper baseline values.
- Mean ~0.094 (exp 001: 0.091, dirichlet: ~0.095).

## Interpretation
Mixing biology + dirichlet does NOT show clear additive gain. The result tracks dirichlet alone, not the sum. This suggests:
- Compositional diversity is the dominant factor that pulled dirichlet > random_uniform.
- DHS biological content is approximately a constant on this harness — doesn't add much beyond what composition + random already provide.
- The ~0.08 ceiling on eval_01 across all strategies tested so far reflects a true noise/signal limit, not just suboptimal sampling.

## Theory updates
- H1 (motifs are core grammar): WEAK evidence so far. DHS sequences (which carry real motifs) don't beat random+composition.
- H5 (dynamic range matters): partially supported. Adding compositional variation helps eval_08 (the hardest eval), suggesting the model benefits from spread along some axis other than pure biology.
- H6 (biology + composition additive): NOT SUPPORTED. The two strategies are not orthogonal in their gains — both seem to access similar variance.
- NEW H7: At 50k, the model is severely under-trained / under-fit. Library design has limited leverage because the trained model isn't expressive enough to capture biological patterns even when present. The dirichlet boost is the model fitting easy GC-content → activity trends, not learning real grammar.

If H7 is right, then bolder strategies might help:
- VERY motif-rich sequences (lots of TF motifs per 200 bp) — give the model big features to fit
- Sequences with EXTREME activity (very strong elements only) — make signal large relative to noise
- Many copies of same elements (low effective N but high precision per element)

## Next experiment idea
EXP 003 — motif-embedded synthetic. Random 200bp backgrounds with 3-8 embedded JASPAR TF motifs per sequence. Tests if heavy motif content (without genomic context) beats DHS. If YES, motifs > genomic context. If NO, the harness can't latch on to motifs at 50k.
