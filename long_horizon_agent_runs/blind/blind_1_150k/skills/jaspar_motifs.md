# Skill: JASPAR motifs for MPRA library design

## Source
JASPAR 2024 CORE vertebrate non-redundant PFMs:
```
https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
```
~1 MB, 4395 lines. Stored at `data/jaspar2024_core_vert_pfms.txt`.
Contains **879** motifs (after parsing); length range 5–33 bp, mean 10.1.

## File format
Each motif is 5 lines: a `>MA<id>.<v> <name>` header, then 4 rows for
A, C, G, T in that order:
```
>MA0004.1   Arnt
A  [   4  19   0   0   0   0 ]
C  [  16   0  20   0   0   0 ]
G  [   0   1   0  20   0  20 ]
T  [   0   0   0   0  20   0 ]
```
The numbers are observed base counts at each position.

## Parser
See `libraries/005_ccre_motif_embedded/generate.py::load_jaspar_consensus`.
Returns consensus DNA strings (argmax base per position). The same parser
is in `libraries/006_ccre_motif_90_10/generate.py`.

If you want probabilistic sampling instead of consensus, normalise each
column to a probability distribution and sample per position.

## Embedding strategy that works (exp 005)
- Random ACGT background, 200 bp.
- 2 motifs per sequence, both placed at random positions (constrained to
  fit in 200 bp).
- Each motif chosen uniformly from the 879 JASPAR motifs.

## Empirical effect
Replacing the uniform-random half of a 75k cCRE + 75k random hybrid with
motif-embedded random gave **+0.003 mean_r** (0.880 → 0.883) across the
14 evals. Tiny but consistent.

The eval_08 super-additivity effect (0.916 in hybrid vs 0.892 cCRE alone)
was *identical* across uniform-random and motif-embedded hybrid sources,
suggesting that effect is about source diversity per se rather than motif
content specifically.

## Ideas to try (untested)
- Probabilistic sampling from PFM instead of consensus.
- More motifs per sequence (5+, denser).
- One motif per sequence (sparser).
- Motif-pair patterns (e.g., two cooperating TF sites spaced 10–30 bp apart).
- Motifs sampled from a curated subset (e.g., only basic helix-loop-helix
  TFs, or only TFs expressed in K562/HepG2/SK-N-SH).
- Strand: half the embedded motifs in reverse complement.
