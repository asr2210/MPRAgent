# 004_dhs_plus_motif_synth — notes

**Design**: 25k DHS uniform + 25k motif-planted synth (3-8 JASPAR motifs each, random backbone).

**Result**: eval_01 = 0.5277, mean = 0.4733.

**Comparison**:
| eval | random | motif | dhs | hybrid |
|---|---|---|---|---|
| 01 | 0.470 | 0.505 | 0.660 | **0.528** |
| 03 | 0.469 | 0.540 | 0.673 | 0.557 |
| 04 | 0.395 | 0.177 | 0.493 | 0.274 |
| 07 | 0.519 | 0.717 | 0.761 | 0.689 |
| 08 | 0.158 | 0.093 | 0.116 | 0.106 |
| 13 | 0.502 | 0.701 | 0.752 | 0.671 |

**The hybrid is WORSE than pure DHS on EVERY eval** (sometimes substantially: eval_04 drops 0.22, eval_01 drops 0.13). H1 (multi-objective spanning) is contradicted.

**Why**: motif-planted synthetic is not "neutral filler" — it actively introduces spurious sequence-to-activity associations that the model learns. These associations don't generalize to real DHS-like test sequences because real motifs are surrounded by real chromatin context (cooperative/competitive TF interactions, modulators, spacing constraints) that planted motifs lack. The model overfits motif-presence signals from the synthetic half, hurting the DHS-derived signal too.

**Contrast with `dhs_synth` baseline (50% DHS + 50% random)** which got 0.7174 — close to pure DHS 0.7089. Random is genuinely neutral filler because it doesn't introduce *any* systematic associations; the model learns nothing from random but also isn't misled. Motif-planted synthetic is *worse than random* as filler precisely because it introduces patterns.

**Theory update (H1 → H2)**:
> Library utility comes from sequences with **biologically plausible joint statistics** spanning the activity range. The bad news: you can't engineer 'pseudo-DHS' sequences by planting motifs into random backbones — the model overfits the planted-motif statistics. The good news: random sequences are accidentally neutral and can be added safely as filler.

**Implications for next experiments**:
- Stop trying to augment with engineered sequences. Need ALL sequences to come from biologically authentic sources.
- The 16 NMF topics cover a wide regulatory program space. Can I get extra value by stratifying across these or by *excluding* the test-cell-aligned topics to test pure generalization?
- The 'extra activity' to fill the low-activity regime should come from biologically real but low-activity regions (deserts, heterochromatin) — not synthetic.

**Next**: test generalization directly. Build a DHS library that EXCLUDES topics aligned with the labeling cell types (K562 ~ Myeloid, HepG2 ~ Digestive/Hepatic, SK-N-SH ~ Neural). If eval_01 stays near 0.66, the library generalizes — labeling-cell-type-specific data is not required, just broad regulatory programs. If it drops a lot, matching matters.
