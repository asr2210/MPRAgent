# 019_ultra_repeat_rich

## Design
50k tiles from chr1/17/19/22 with extreme repeat-rich filter:
- max single-nuc > 0.55
- OR 2-mer entropy < 2.5
- OR homopolymer ≥ 8
Acceptance rate: 12% (50k of 415k attempts). seed=0.

## Result
- eval_01 = **0.4698** (vs E2 0.4992, E18 0.4937; Δ -0.030 vs E2)
- mean over 14 ≈ 0.466

## Interpretation
Confirms monotonic narrowing penalty. The complexity-axis pattern:
- All random (E2): 0.499
- Repeat-rich 70%: 0.494
- Ultra-repeat 12%: 0.470
- High-complexity 30% (E17): 0.407

Peak is at "almost all random" — anything narrower hurts. The
informative subset isn't a small subset; it's the BROAD complexity
range. Even my "informativeness ranking" of E18 vs E17 was over-
interpreted: yes high-complexity tiles contribute LESS per slot,
but every additional slot of any reasonable complexity helps.

## Theory v8 confirmed

**The information content of random human DNA at 50k is roughly the
broad complexity distribution itself.** No single complexity tier is
the "main" signal — the eval depends on diversity across the tier.

The 0.50 ceiling at 50k = the entropy of the broad complexity prior.
Any filter (in either direction) reduces effective sample diversity
without proportional gain.

## What I've ruled out

After 19 experiments, the following do NOT break 0.50:
- DHS variants (stratified, uniform, top-signal): 0.44-0.47
- Cross-species (mouse alone or mixed): 0.45-0.48
- Mixed sources (DHS+genome, mouse+human, top-DHS+random): 0.48-0.50
- Element-class balance (cCRE PLS/pELS/dELS/other): 0.40
- Promoter-only / TSS-centered: 0.37
- Synthetic (uniform, Markov, motif-implant): 0.27-0.35
- Filtered random (high-complexity OR ultra-repeat): 0.41-0.47

What's left untested:
- **Gene-density weighted across all chroms** — most natural broadening
- **GC-balanced random** — tests GC distribution match
- **Small-fraction mixtures** (10% something + 90% E2)
- **Synth + real hybrids** (could give regularization)

## Plan for E20
**Gene-density-weighted broad sampling**: compute gene count per 1Mb
bin across all 24 chroms, sample tiles weighted by local gene
density. Tests if continuous gene-density weighting beats hard
chromosome filter (E2 vs E6 gap was 0.499 vs 0.482).

Predicted: 0.49-0.50. Probably saturated.

If >0.50: minor win, push gene-richness more.
If ≈0.50: confirms natural ceiling.
If <0.49: diversity gradient hurts somehow.
