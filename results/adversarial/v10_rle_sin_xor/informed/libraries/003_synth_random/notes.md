# Experiment 003 — Synthetic random (diagnostic)

## Hypothesis
If synth_random gives SK-N-SH > 0, then the issue in exp 001/002 was DHS-specific.
If synth_random also gives SK-N-SH ≈ 0, the SK-N-SH issue is pipeline-wide.

## Result
eval_01 = 0.5210. K562=0.99, HepG2=0.56, SK-N-SH=0.

This matches the `random_uniform` baseline in strategies.md (0.5202) — strong evidence
that strategies.md is the relevant baseline reference, NOT instructions.md Table 1.

## Verdict
SK-N-SH≈0 is a SYSTEMATIC pipeline characteristic, not DHS-specific.
K562≈0.99 for random sequences indicates K562 oracle is trivially learnable
for any library. The headroom for mean_r improvement is in HepG2 and SK-N-SH.

## Insight
- My pipeline matches strategies.md baselines, not instructions.md Table 1
- random_uniform=0.5202 is the strongest known simple baseline → my target to beat
