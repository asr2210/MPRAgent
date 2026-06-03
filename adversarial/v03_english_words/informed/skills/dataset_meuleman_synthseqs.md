# Skill: Meuleman synthseqs DHS dataset

## What it is
A pre-extracted, pre-curated subset of the Meuleman 2020 DHS Index, packaged for
training a 16-class classifier of NMF components (cellular contexts). Available at:
- https://www.meuleman.org/train_all_classifier_light.csv.gz (160k seqs, 15MB)
- https://www.meuleman.org/validation_all_classifier_light.csv.gz (16k seqs, chr2)
- https://www.meuleman.org/test_all_classifier_light.csv.gz (16k seqs, chr1)

## Format
TSV (gzipped, despite `.csv` extension). Columns:
- `seqname`, `start`, `end`, `DHS_width`, `summit`
- `total_signal` (DHS signal magnitude)
- `numsamples` (in how many of 733 biosamples this DHS is called)
- `raw_sequence` (200bp, ACGT only, summit-centered, always exactly 200 chars)
- `component` (1..16, the dominant NMF cellular-context component)
- `proportion` (NMF loading proportion in that dominant component)

Training set: 10k per component, sourced from chr3..chrY.

## Critical bias (caused 001 to fail)
**This is NOT a representative DHS sample.** It is curated for cell-type
**specificity**:
- Mean `proportion` = 0.87, median = 0.94 (highly specific to one component)
- Median `numsamples` = 6 of 733 (rare elements)
- 50% of sequences active in ≤6 biosamples

So sampling uniformly (or stratifying by component) from this pool produces a
library *heavily biased toward rare, cell-type-specific elements that are likely
INACTIVE in K562/HepG2/SK-N-SH*. Most sequences will have ≈zero MPRA signal in our
three measured cell types, which gives the model nothing to learn.

Evidence: experiment 001 (NMF-stratified × signal-quartile-stratified) got
eval_01 = 0.388 vs dhs_topic baseline of 0.7232. SK-N-SH correlation was ~0.06
across all evals.

## How to use it sensibly
- **Bias toward high `numsamples`**: sequences detected in many biosamples are more
  likely to be active in K562/HepG2/SK-N-SH (which gives the model real signal).
- **Bias toward high `total_signal`**: signal-strong elements have measurable
  activity → trainable labels.
- **Be careful with NMF-stratified sampling**: this enforces equal coverage of
  cellular contexts, several of which may have ZERO overlap with our three
  labeling cell types.
- **For a representative sample of all DHSs**, use the full DHS Index
  (~3.6M elements, requires hg38 FASTA to extract sequences) — not this subset.

## Loading snippet
```python
import pandas as pd
df = pd.read_csv("data/train_all_classifier_light.csv.gz", sep="\t")
# Each df["raw_sequence"] is exactly 200bp, ACGT.
```
