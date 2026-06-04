# 003_motif_packed

## Design
50,000 synthetic 200bp sequences. Each contains 4-8 strong consensus TF
motifs from a curated list of 17 (AP-1, CREB, NF-kB, SP1, TATA, E-box,
NF-Y, GATA, ETS, C/EBP, HNF1/4, KLF, OCT, FOXA, etc.) placed at random
non-overlapping positions in a random ACGT background.

## Hypothesis
If MPRA activity is driven by motif content and the v14 model can learn
even simple motif→activity rules in brief training, this library should
give the model the cleanest possible signal: every motif appears in many
sequences, model can learn each motif's coefficient.

## Result
eval_01 mean_r = -0.0037. Still in the v14 noise floor (-0.004 to +0.003).

## Interpretation
Motif content alone — even tightly packed and recurring — does not get
v14 above its noise floor. Either the model is not learning motif→activity
mappings at all, or the eval sets do not value motif presence.

What this rules out:
- Simple motif counting is not the v14 model's mechanism.
- Single-class libraries (all the same kind of sequence) do not work,
  regardless of whether the class is biological (cCRE), motif-rich, or
  random.

Next idea: bimodal library — half motif-rich, half motif-poor — to give
the model strong label variance to learn from.
