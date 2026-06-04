# Skill: Loading hg38 and the Meuleman DHS index

## Files

- `data/hg38.fa` (~3.3 GB) — full hg38 reference. Download:
  ```
  curl -L -o data/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip data/hg38.fa.gz
  ```
- `data/dhs_index.txt.gz` (~90 MB) — Meuleman 2020 DHS index, 3,591,898
  elements with NMF-component labels (16 topics). Download:
  ```
  curl -L -o data/dhs_index.txt.gz \
    https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
  ```
- `data/hg38.pkl` — pickled `{chrom: uppercase_str}` cache built on first
  call to `load_genome()`. Loads in ~30 s vs ~3 min for re-parsing.

## API (skills/genome_utils.py)

```python
from skills.genome_utils import load_genome, load_dhs, fetch, extract_200bp_around_summit

genome = load_genome()          # {chrom -> 'ACGT...'} ; cached pickle
dhs = load_dhs()                # pandas.DataFrame
seq = fetch(genome, "chr1", 16170 - 100, 16170 + 100)  # 200 bp
seqs = extract_200bp_around_summit(genome, dhs.sample(50_000), 200)
```

## DHS index columns

`seqname, start, end, identifier, mean_signal, numsamples, summit,
core_start, core_end, component`

- `summit`: high-confidence center; safer than midpoint for windowing.
- `numsamples`: how many biosamples this DHS is open in (1 = ultra-cell-
  type-specific, 733 = tissue-invariant).
- `component`: one of 16 NMF topic labels (e.g. "Tissue invariant",
  "Placental / trophoblast", "Lymphoid", "Cardiac", ...). Topic loadings
  themselves are not in the index; `component` is the argmax topic.
- `mean_signal`: average normalized DNase signal across positive samples.

## Gotchas

- Some DHS rows extract sequences with `N`s near assembly gaps — drop
  these (`extract_200bp_around_summit` returns empty string for those).
- 200 bp centered on summit may extend past chromosome ends for ~tens of
  rows; same filtering catches it.
- `dhs_topic`-style sampling uses the row's `component` argmax — Meuleman
  also publishes full NMF loadings as a separate file. For most weighting
  schemes the argmax is sufficient; if soft topic weighting is needed,
  pull the loadings file separately.

## Sanity numbers (first time using this)

- 3.59 M total DHS rows
- ~3.5 M survive 200 bp extraction with no Ns
- Largest components: "Tissue invariant" (~1.0 M), "Stromal A" (~0.4 M),
  ... ; smallest: "Erythroid" (~0.04 M) — important when stratifying.
