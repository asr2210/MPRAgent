# 001_ccre_random

## Design
Uniform random sampling of ENCODE cCREs (V3, GRCh38), 1.06M elements.
For each sampled cCRE, extract 200bp centered on midpoint. Reject any window
crossing chromosome bounds or containing N. 50,000 sequences total.

## Hypothesis
Was meant to be a clean reproduction of the dhs_random baseline (eval_01 ≈
0.7089 per instructions.md Table 1). cCREs are an ENCODE-curated subset of DHS,
so I expected r in the same ballpark.

## Result
mean_r on eval_01 = **0.0016**. All eval sets near zero (range -0.0019 to
0.0054). Total time 11.6s, n_seeds=1.

## Interpretation — major surprise
The "informed" baseline table in instructions.md predicts ~0.71 for this kind
of library. The actual v14 baselines in strategies.md show all simple synthetic
strategies getting ~0. My biologically-curated cCRE library also gets ~0.

This means either:
- (a) the instructions.md baseline table does not actually apply to this v14
  prepare.py (the informed instructions may be intentionally misleading), or
- (b) cCREs differ from DHS in a way that matters, or
- (c) something pathological about my sequence extraction.

Sequence inspection looks clean (52% AT/48% GC, complex sequences, all 200bp
ACGT). Most likely (a) — v14's evaluator behaves very differently from what
instructions.md claims.

n_seeds=1 and total time 11.6s suggest v14 uses an extremely fast / minimal
training process. With a single seed and tiny compute budget, the model may
only be able to learn libraries with very strong, very simple signal.

## Implication for theory
My v0 theory ("diverse biological regulatory sequences will generalize") is
not predictive of v14's behavior. Diverse biology gave no signal. I need to
think about what kind of library can teach a quickly-trained model anything
useful at all in v14.
