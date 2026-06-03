# Experiment 001 — dhs_signal_specificity

## Design
Sample 50,000 Meuleman 2020 DHS elements with weight ∝ `mean_signal /
sqrt(numsamples)`, extract 200 bp centered on the summit, no filtering on
component.

## Hypothesis
Should land within ±0.02 of the published `dhs_topic` baseline
(`eval_01 = 0.7232`), since both target "elements with strong cell-type-
specific accessibility signal".

## Result
`eval_01 = 0.5600` — **0.16 lower** than the dhs_topic baseline, and
even ~0.15 below `dhs_random` (0.7089). Pattern across evals is flat
(~0.52–0.60) except `eval_08 = 0.1414` (matches synth-collapse signature
seen in dhs strategies generally).

## What this tells me
The signal/numsamples weighting biases way too hard toward elements with
`numsamples == 1`. The 25th percentile of numsamples is literally 1, and
my weight gives those the highest probability — so the bulk of my draws
are single-biosample peaks, which are exactly the noisiest, least-
reproducible DHS calls. The model can't learn from those.

**Updated understanding:** "Cell-type-specific" in the dhs_topic baseline
description refers to *topic-loaded* sequences in the NMF decomposition —
which concentrates on elements that have a strong, *interpretable* program
membership, not just elements that happen to be in one biosample. Topic
loadings are smoothed across correlated biosamples; raw `1/numsamples`
isn't. So `1/numsamples` is a bad proxy.

## Next
Calibrate by running uniform-random DHS sampling — that should land near
the `dhs_random` baseline (0.7089). If it does, the harness and my data
extraction are fine and only the weighting was wrong. If it doesn't, the
data source or window definition is wrong and I need to investigate.

## Numbers
mean_r averaged across 14 evals: 0.510
eval_01: 0.5600
eval_08: 0.1414
time_s: 16.4
