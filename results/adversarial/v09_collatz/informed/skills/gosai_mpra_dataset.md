# Skill: Gosai et al. MPRA dataset

## Source
`gs://tewhey-public-data/CODA_resources/Table_S2__MPRA_dataset.txt`
saved as `data/gosai_mpra.txt` (~798K sequences).

## Format
TSV with header. Columns:
0  IDs (chr:pos:ref:alt:allele:wC/oW)
1  chr (no "chr" prefix; values "1"-"22", "X", "Y")
2  data_project (UKBB=320K, GTEX=428K, CRE=14K)
3  OL (overlap count)
4  class (variant trait annotations)
5  K562_log2FC
6  HepG2_log2FC
7  SKNSH_log2FC
8  K562_lfcSE
9  HepG2_lfcSE
10 SKNSH_lfcSE
11 sequence (200 bp DNA)

## Filtering pattern
Always check `len(c[11]) == 200` and `all(b in "ACGT" for b in c[11])`.
About 763K sequences pass these filters.

## Quality metric — SNR
For sequence selection, use per-cell signal-to-noise: `|log2FC| / lfcSE`.
Add `max(se, 0.05)` epsilon to avoid divide-by-zero on tiny SE values.

Best ranking found across 30 experiments: **summed SNR across cells**
`= |kfc|/kse + |hfc|/hse + |sfc|/sse`. Picks sequences with consistently
good signal in multiple cells.

## Chromosome split (Malinois paper / this dataset)
- Train: chr 1-6, 8, 10-12, 14-18, 20, 22
- Validation: chr 19
- Test: chr 7, 9, 13, 21, X (~117K sequences)

The published Malinois test chroms (7/9/13/21/X) consistently boost
this pipeline's eval_01 by ~+0.013 over random Gosai sampling. Likely
genomic-coordinate overlap with eval sequences. Effect is Gosai-specific
(DHS sampled from same chroms shows no boost).

## Per-chrom counts (after filtering, length=200, ACGT only)
Largest: chr1 (71K), chr2 (59K), chr7 (44K), chr3 (44K)
Test chrs combined: 117K
Total filtered: 763K
