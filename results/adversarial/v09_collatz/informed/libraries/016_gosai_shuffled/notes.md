# E016 — Gosai shuffled (base-level scramble per sequence)

Take 50K Gosai sequences, shuffle bases within each sequence. Preserves
mononucleotide composition (GC%, A%, etc.) but destroys all dinucleotide
structure, motifs, and positional grammar.

## Result
eval_01 = 0.3052. Per-cell K562=0.140, HepG2=0.188, SKNSH=0.588.

## Interpretation — major finding
This sits BETWEEN E001 random uniform (0.244) and E004 Gosai random (0.323).

Decomposition of Gosai's +0.079 lift over random uniform:
- ~+0.061 from base composition alone (matching natural GC)
- ~+0.018 from real motif grammar / sequence structure

77% of the Gosai-vs-random gap is just **base composition**. Real
regulatory grammar contributes ≤+0.02 to my pipeline's eval_01. This
is consistent with the K562/HepG2 ceiling story: the pipeline can't
efficiently learn from motif-level signal, only from gross composition.

## Implication
Library-design effort focused on motif content or sequence-grammar
sophistication is mostly wasted. The big lever is matching natural
base-frequency distribution. The ceiling at 0.34 likely reflects the
small remaining grammar component.
