# Skill — Meuleman DHS curated training set

A pre-extracted 200bp DHS sequence dataset, perfect for quick DHS-based experiments WITHOUT needing hg38 + DHS BED + sequence extraction.

## Source

- Training: https://www.meuleman.org/train_all_classifier_light.csv.gz (15 MB, 160k seqs)
- Validation: https://www.meuleman.org/validation_all_classifier_light.csv.gz (1.5 MB, chr2)

## Schema (tab-separated csv inside .gz)

Columns: `seqname start end DHS_width summit total_signal numsamples raw_sequence component proportion`

- `raw_sequence`: 200bp uppercase ACGT (verified — all 200bp, no Ns in training set)
- `component`: 1-16 (NMF component majority label)
- `proportion`: how confidently the majority component dominates (0.35-1.0, mean 0.87)
- `total_signal`: aggregate accessibility signal across all 733 biosamples
- `numsamples`: number of biosamples in which the DHS was called
- Training set: chr3-chrY. Validation: chr2. Test (not used): chr1.
- Training is balanced: exactly 10,000 per component.

## Loading

```python
import pandas as pd
df = pd.read_csv("data/train_all_classifier_light.csv.gz", sep="\t")
# df.raw_sequence — 200bp ACGT strings, ready to write
# df.component   — integer 1..16
```

## Sampling recipes

**Stratified (equal per component)**: `df.groupby("component").sample(n=N//16, random_state=0).raw_sequence`

**Topic-weighted (proportional, but already balanced so this collapses to uniform within component)**: weight by proportion within component, then take all.

**Signal-weighted**: `df.sample(n=N, weights=df.total_signal, random_state=0).raw_sequence`

**By num biosamples**: `df.sample(n=N, weights=df.numsamples, random_state=0).raw_sequence` — biases toward broadly accessible (housekeeping-like) DHSs.

## Caveats / gotchas

- 160k is small relative to the full 3.6M DHS index. The full index gives more diversity (e.g., for topic-proportional sampling).
- Training set excludes chr1 and chr2, so spatial bias from chromosome composition exists.
- The "components" are Meuleman NMF over 733 biosamples; not the same as ChromHMM states or SEI states.
- Some downstream evaluators may have eval sets that overlap with these training sequences (unknown but possible).

## Empirical result (exp 001 in this run)

Sampling 50k stratified (3125/component) from training set → eval_01 = 0.0739 in this harness. Random_uniform baseline = 0.0711. So in this harness, DHS-stratified barely beats random — the realistic ceiling here looks much lower than instructions.md Table 1 suggests.
