# 006_dhs_nolabel_plus_random — notes

**Design**: 25k DHS from 12 non-label-aligned topics + 25k uniform random. Tests whether stacking exp 005's exclusion with `dhs_synth` baseline's random filler helps.

**Result**: eval_01 = 0.5349 — drops 0.14 vs exp 005 (0.6752). On every eval, exp 006 is worse than exp 005.

| eval | exp 005 (50k DHS-NL) | exp 006 (25k DHS-NL + 25k rand) | Δ |
|---|---|---|---|
| 01 | 0.6752 | 0.5349 | -0.140 |
| 04 | 0.5374 | 0.3201 | -0.217 |
| 07 | 0.7597 | 0.6708 | -0.089 |
| 13 | 0.7509 | 0.6516 | -0.099 |

**Stacking failed catastrophically.** In MY pipeline, replacing half the DHS budget with random ACTIVELY HURTS — opposite to what the `dhs_synth` baseline (0.7174 ≈ dhs_random 0.7089) implied.

**Why my pipeline differs from baseline**: the baseline's `dhs_synth` may use a different labeling scheme (e.g. partial oracle labels for the synth half, or oracle labels for everything). My pipeline gives noisy MPRA labels to all 50k sequences. Random sequences carry near-zero per-sequence information; the model wastes capacity learning their noise; the DHS half's effective grammar learning is diluted.

**Theory update H2 → H3**:
> *In this prepare.py pipeline*, training-data utility scales with informative-label signal density. Mixing in low-information sequences (random) dilutes the effective dataset. **Pure-source libraries beat mixtures.** Exp 004 (DHS + motif-planted) and exp 006 (DHS + random) both showed this; the baseline table's apparent support for mixtures (`dhs_synth` ≈ `dhs_random`) is likely an artefact of a different label source for the synth half in that experiment.

**Implication**: stop mixing. Pursue PURE-source library variants. The best library so far is exp 005 (DHS from non-label-aligned topics). Push that direction by improving DHS sampling quality:
- Quality weighting (mean_signal × numsamples)
- Topic-balanced sampling within non-label-aligned pool
- Motif-density selection within DHS

**Next**: pick one direction. Most promising: **`numsamples`-weighted DHS sampling from non-label-aligned topics**. Tests whether biologically robust (broadly accessible) DHSs are better training material than narrow/single-biosample DHSs.
