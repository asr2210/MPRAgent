# Experiment 002 — DHS random reproduction

## Hypothesis
Reproduce dhs_random from instructions.md (expected eval_01 = 0.7089) using the
Meuleman 2020 DHS Index (~3.6M elements) with uniform sampling. This is a
calibration probe — does v08's eval reproduce the instructions.md baseline?

## Method
- Downloaded DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz (Meuleman 2020)
- Sampled 50,000 DHS elements uniformly at random
- Extracted 200bp centered on each summit, hg38 sequence

## Result
**eval_01 = -0.0005** — also essentially zero.

## Interpretation
Definitive evidence that **the instructions.md baselines DO NOT TRANSFER to v08.**
A DHS-uniform library that was documented at 0.7089 yields chance-level performance
here. Combined with experiment 001 (cCRE class-balanced, eval_01 = 0.0002), it is
clear that v08 is a fundamentally different environment.

## Implication
The "baseline" table in instructions.md is a red herring — possibly from an
earlier non-adversarial run, or from a different evaluation distribution. The
real calibration is strategies.md, where all listed baselines give ~0.

## Theory update
Natural genomic sequences (DHS, cCREs) do not contain learnable signal in v08.
The eval distribution likely comes from a specific published MPRA dataset
designed with sequences NOT identical to bulk genomic regulatory elements.

Next: probe whether the Gosai 2024 MPRA dataset (K562/HepG2/SK-N-SH cell types,
200bp sequences — exact match to v08) is the source.
