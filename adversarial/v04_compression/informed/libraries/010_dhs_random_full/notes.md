# 010_dhs_random_full

## Design
50k DHS sites sampled uniformly at random from the full Meuleman DHS
index (~3.6M elements on chr1–22+X+Y). For each: 200bp window centered
on summit. seed=0.

This is the `dhs_random` strategy: no NMF stratification, just uniform
sampling of DHS sites.

## Result
- eval_01 = **0.4699** (vs E2 random genome 0.4992; E3 dhs_strat 0.4378)
- mean over 14 ≈ 0.476

Clean middle ground:
- DHS-stratified (E3, 0.4378) < DHS-random (E10, 0.4699) < random genome (E2, 0.4992)
- Δ (E3 → E10): +0.032 — removing NMF-stratification helps.
- Δ (E10 → E2): +0.029 — removing DHS-restriction helps.

Both axes have ~0.03 r-point penalties. Both narrowings hurt.

## Theory update v4
Two distinct narrowness penalties:
1. **Element-type narrowing** (using DHS-only vs genome-broad): -0.03 r.
2. **NMF stratification** (forcing equal coverage across cell-context
   programs): -0.03 r. Stratification over-weights rare/specialized
   programs which carry signal less aligned with the broad eval distribution.

Composite: an aggressively-curated library (DHS + NMF-stratified) is
0.06 r below broad random genome. Each curation step costs ~0.03 r.

## Implications
- For design: avoid stratification on functional axes; let the
  natural distribution dominate.
- For theory: the eval distribution is closest to broad random hg38,
  with mild gene-richness preference (chr1/17/19/22 > all chroms).
- The published baseline's claim that DHS dominates synth doesn't
  hold in this pipeline; instead, "real DNA dominates synth" is the
  big effect, and curation within real DNA mostly hurts.

## Confidence check
This is consistent with E9 (promoter-only = 0.37) which is the
extreme narrowing case. The progression:
- random hg38 broad (gene-rich):  0.499
- random hg38 broad (all chroms): 0.482
- DHS random (uniform):           0.470
- DHS stratified by NMF:          0.438
- Promoter TSS only:              0.370
- 4-mer Markov synth:             0.268
- uniform random synth:           0.307

Strong monotonic pattern. Narrowness is bad.

## What this DOESN'T tell me
This rules out a class of designs (anything that narrows). To go
ABOVE 0.50, I need designs that don't narrow but somehow *add
information* — through diversity expansion (cross-species, repeat
classes, mixed sources) or activity stratification using an oracle.
