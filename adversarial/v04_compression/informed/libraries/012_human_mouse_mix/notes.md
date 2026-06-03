# 012_human_mouse_mix

## Design
25k random tiles from human chr1/17/19/22 + 25k from mouse mm10
chr1/11/17/19. Shuffled. seed=0.

## Result
- eval_01 = **0.4839** (vs E2 0.4992 human, E11 0.4485 mouse)
- mean over 14 ≈ 0.483
- Linear mix expectation = (0.4992+0.4485)/2 = 0.4739
- Mix beats linear by +0.010 (small diversity gain)
- Mix loses to pure human by -0.015 (dilution dominates)

Across every eval: pure_human > mix > pure_mouse. The dilution
penalty (~5 r-points cost from mouse) exceeds the diversity benefit
(~1 r-point gain from cross-species).

## Interpretation
The model uses mouse DNA effectively (E11 = 0.45, far above random),
but mixing it with human in equal proportion costs more than it adds.
The eval distribution is human-biased; sequences that are NOT human
provide weaker signal per slot.

This is the second confirmation of "mixing in non-pure-human sources
slightly hurts":
- E4: 25k human + 25k DHS-stratified = 0.4926 (vs pure human 0.4992)
- E12: 25k human + 25k mouse = 0.4839 (vs pure human 0.4992)

Both ~0.01-0.015 below pure human. Pure broad-genome-human is a
local optimum.

## Theory update v5.1
"Mixing strategies" with 50/50 splits don't break the human ceiling.
Two possible escapes:
(a) Smaller fraction of admixture (10-20% non-human) — limit dilution.
(b) Mix in something OTHER than DHS or mouse — must be incrementally
    informative, not just real DNA from a different distribution.

The "cross-species adds info" hypothesis is partially true (+0.01)
but cross-species also dilutes (-0.025 per unit of mouse). At 50/50
the dilution wins.

## What to try next
**Selective injection:** 40k human random + 10k high-activity DHS
sites (top mean_signal across 733 biosamples). The mean_signal metric
ranks DHS sites by universal-cell-type activity — these are the
sequences most likely to generalize across UNSEEN cell types.

Prediction: if the eval set is testing "transfer to new cell types",
high-mean-signal DHS sites should be the best supplement. 20%
replacement is small enough to avoid the E4/E12 dilution disaster.

Outcomes:
- ≥0.51: high-signal injection adds genuine cell-transfer info.
- ~0.50: neutral, info gain ≈ dilution.
- <0.49: even small DHS injection hurts (decisive: stop adding non-
  random-genome things to a random-genome backbone).
