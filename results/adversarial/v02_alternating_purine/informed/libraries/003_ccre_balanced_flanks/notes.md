# Experiment 003 — cCRE 5-class balanced + paired flanks

## Design
- 25k positives: 5,000 each from PLS, pELS, dELS, CTCF-only, DNase-H3K4me3
- 25k paired flanks (±1500-3000bp, no cCRE overlap)
- Tests: does class-balanced beat blind's uniform cCRE+flanks (0.166)?

## Result — eval_01 = 0.1445 (WORSE than blind 010 = 0.166)
- K562_r = -0.009 / SKNSH_r = 0.451
- SKNSH is excellent (highest of any exp so far) but K562/HepG2 dropped

## Key finding
Forcing class balance HURTS. Natural cCRE distribution is 74% dELS, 16% pELS,
4% PLS, 3% CTCF-only, 2% DNH3. My 20%/20%/20%/20%/20% over-represents rare
classes and under-represents the dominant dELS signal.

The model presumably learns dELS-grammar from many examples (74k → 5k = 15x
fewer than natural). Rare classes may saturate quickly (CTCF doesn't need 5k
examples to learn).

## Theory update (Theory v3)
**Equal naming ≠ equal information.** Library composition should track
*informational return per class*, not nominal balance. Blind 013's asymmetric
recipe (keep dELS natural, BOOST rare from <1% to 10%) works because it
preserves dELS signal while ensuring rare classes are merely visible.

## Tower of evidence
- Exp 001: DHS topic = 0.129 (source: DHS, recipe: weighted)
- Exp 002: DHS topic + flanks = 0.133 (+flanks ≈ +0.004)
- Exp 003: cCRE 5-balanced + flanks = 0.144 (cCRE > DHS but balance hurts)
- Blind 010: cCRE uniform + flanks = 0.166
- Blind 013: cCRE asymm + flanks = 0.173

## Open question
Blind 013 noted eval_13 dropped because PLS was starved (only ~570 in 15k
uniform). Adding PLS boost may help eval_13 without hurting other evals.
