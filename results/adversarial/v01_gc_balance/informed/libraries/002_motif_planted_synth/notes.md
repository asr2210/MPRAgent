# 002_motif_planted_synth — notes

**Hypothesis tested (H0a)**: regulatory grammar density drives library utility. If planting JASPAR TF motifs into otherwise-random sequences substantially improves over uniform random, motif identity is what the model picks up.

**Design**: 50k × 200bp. Each sequence: random backbone + 3-8 motifs from JASPAR 2024 vert non-redundant (879 motifs), random positions, random strand, instances sampled from each PWM.

**Result vs exp 001 (random uniform)**:
| eval | random | motif-planted | Δ |
|---|---|---|---|
| 01 | 0.4704 | 0.5052 | +0.035 |
| 02 | 0.4702 | 0.5071 | +0.037 |
| 03 | 0.4688 | 0.5395 | +0.071 |
| 04 | 0.3947 | 0.1773 | **-0.217** |
| 07 | 0.5192 | 0.7171 | **+0.198** |
| 08 | 0.1576 | 0.0925 | **-0.065** |
| 09 | 0.3947 | 0.1773 | **-0.217** |
| 13 | 0.5020 | 0.7009 | **+0.199** |

Mean: 0.4433 → 0.4642 (overall +0.02, very small).

**This is a dramatic reshape**, not a uniform lift. Motif planting helped eval_07 and eval_13 (gain +0.20) while badly hurting eval_04, eval_09 (drop -0.22). Interestingly, eval_04 ≡ eval_09 (identical pair); eval_07 and eval_13 are unique but both gained similarly.

**Interpretation**: 
- eval_07, eval_13 evidently reward the model's ability to recognise TF motif presence. They are likely sequences engineered or selected for clear motif content. Motif planting in training makes the model more sensitive to motif occurrence → big gain.
- eval_04, eval_09 (identical pair) presumably evaluate moderate-activity baseline sequences. Adding salient motifs *shifts* the training distribution toward extreme activities, hurting calibration on moderate-activity test sequences. Possibly these eval sets emphasise the *absence* of motifs (insulators, decoys).
- eval_08 (already very hard for everyone) gets even worse — its sequences may be specialised in a way that motif-planted noise actively hurts.

**Theory update (H0 → H0a + new finding)**: 
Library design is **multi-objective by eval set**. A single library cannot simultaneously maximise all evals: there is an apparent **motif-content axis** that helps motif-dependent evals (07, 13) and hurts baseline-calibration evals (04, 09).

**Implication**: A good general library probably needs to span BOTH ends — sequences with strong motifs (e.g. DHS, planted motifs) AND sequences without (random, scrambled, neutral controls). This explains why `dhs_synth` (DHS + 50% random) is roughly tied with `dhs_topic` despite half the sequences being "wasted" on random — the random half is *not* wasted, it's covering the low-motif regime.

**Next**: build DHS extraction pipeline (exp 003) and confirm I can reach ~0.72 on eval_01. Then I can start mixing strategies informed by the multi-objective insight.
