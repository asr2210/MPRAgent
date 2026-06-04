# Experiment 001 — dhs_topic reproduction

## Design
- 50k sequences from DHS pool (Meuleman et al. 2020, 3.59M elements)
- Weighted by sum of NMF topic loadings (16 topics) per element
- Single seed (0)
- Intended as my pipeline anchor against the published dhs_topic = 0.7232

## Result — eval_01 = 0.1293
mean_r ∈ [0.048, 0.169] across 14 evals.
K562_r and HepG2_r are tied to ~0 across all evals — model collapsed.
SKNSH_r ≈ 0.4 carries the signal.

## Major finding
The published baselines do not apply to this harness. Inspecting v02
blind_claude runs (same harness) shows EVERY result has K562_r == HepG2_r
identical to 4 decimals across all 14 evals. The harness has structurally
tied K562/HepG2 channels and the achievable range is ~0.12–0.17, not 0.7+.

## Implication
The 30-experiment search is in a different optimization landscape than the
baseline table suggests. I need to think about what helps SKNSH (and any
K562 differentiation I can get) rather than what would have helped a
stronger Malinois-style model.
