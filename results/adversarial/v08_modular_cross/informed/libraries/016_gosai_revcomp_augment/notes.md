# Experiment 016 — Reverse-complement augmentation

## Method
25K Gosai random (lfcSE<0.5) + 25K revcomps of the same 25K. Tests whether
recycling measurements via RC adds training signal.

## Result
**eval_01 = 0.0102** — below 0.018 plateau. Splitting unique sequences in half
to fit RC pairs costs too much diversity.

Notable: **eval_13 = 0.0101** — highest yet for that eval (others typically
≤0.005). RC augmentation specifically helps eval_13. eval_13 may test
"unnatural" sequences that benefit from the model seeing both strands.

## Interpretation
Sequence diversity > measurement repetition. The oracle gives ~equivalent
labels to revcomps (if it's RC-symmetric), so RC adds no new information per
slot. Better to use those slots for new measurements.

But the eval_13 lift is real and suggests some eval sets specifically reward
"both-strand training" — a discriminative signal.

## Next directions
- Test if GTEX-only random+quality beats mixed (017, parallel to 014 for GTEX)
- Test Malinois-training-chromosome subset (018) — if oracle accuracy depends
  on chr provenance
