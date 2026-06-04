# Experiment 011 — Mode-stratified Dirichlet(0.5)

## Result
eval_01 = **0.1356** (vs dirichlet(0.5) 0.1395 — WORSE by 0.004)

## What I did
Generated 500k Dirichlet(0.5) candidates, categorized by
(dominant base, second-dominant base) — 12 modes (4*3). Sampled
~4166 per mode to guarantee uniform coverage of (dom, sec) modes.

## Interpretation
Forcing uniform coverage across 12 composition modes HURTS slightly.
Natural Dirichlet(0.5) sampling already over-represents the "interior"
of each mode (where dom-sec ratio is moderate), which matches the
evaluator preference. Forcing uniform mode coverage adds extra weight
to corner-like compositions, similar to alpha<0.5 — and we know that's
worse (alpha=0.3 → 0.1344).

## Plateau confirmation
Three different alpha-related experiments now all hover at 0.1344-0.1377:
- alpha=0.3 → 0.1344
- mode-stratified → 0.1356
- sobol uniform → 0.1371
- mixed alpha 0.3-2.0 → 0.1377
- alpha=0.5 → 0.1395 (still best)
- alpha=1.0 (=sobol) → 0.1371

The Dirichlet(0.5) sweet spot is very robust. Single-seed noise ~0.005.

## Next
Try a fundamentally different axis: cell-type targeted synthseqs.
K562 head correlation is ~0.04 (much lower than HepG2/SK-N-SH at
~0.17/0.20). Mean improvement might be best driven by lifting K562.
SynthSeqs component 15 (Myeloid/erythroid) is K562-relevant.

**Exp 012**: Topic-uniform synthseqs (3125 per NMF component × 16
components = 50k). vs exp 001 which used natural component frequency
(over-weighted embryonic, under-weighted myeloid). Tests whether
uniform topic balance helps K562 or just shifts the loss elsewhere.
