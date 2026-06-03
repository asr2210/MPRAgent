# Experiment 001 — DHS signal-weighted, summit-centered

## Hypothesis
DHS regions weighted by `mean_signal × log(numsamples+1)` and centered on summit
will roughly reproduce the dhs_topic baseline (~0.72 eval_01).

## Setup
- 50,000 sequences, 200bp, single seed
- Weighted draw from 3.59M DHS elements (Meuleman 2020)
- Standard primary chromosomes only

## Result
- eval_01 mean_r = **0.493** — far below dhs_topic baseline (0.7232)
- eval_08 = 0.4365 — also below baselines
- Striking per-cell-type pattern: K562=0.92, HepG2=0.56, **SK-N-SH≈0**

## Interpretation
Mean-signal weighting heavily concentrates on a small number of very high-signal
elements. The signal distribution is extremely long-tailed (median 0.41, max 427).
Even with mild numsamples boost, the top 20% of weight density falls on the top 50k
elements — so this is effectively a top-k selection of bright regions.

Bright DHS regions appear to be dominated by:
- Strong K562 elements (great for K562 prediction)
- Common/tissue-invariant elements
- Little SK-N-SH-specific signal

The model thus learns K562 well, HepG2 partially, and SK-N-SH not at all.

## Key Learning
**For cross-cell-type generalization, signal-weighted sampling is harmful.**
It biases toward "common/bright" regions that don't teach cell-type-specific
regulatory grammar. Cell-type DIVERSITY in the library matters more than
average signal strength.

## Next Experiment
Reproduce dhs_random (uniform DHS sampling) as a properly-calibrated baseline.
If that matches the reported 0.7089, my pipeline is sound and I can iterate from
there with confidence.
