# Skill: Available reference data

All under `/data/users/arao/sandbox_runs/blind_2_150k/data/` — symlinked into
`./data/` in this run. **Use `data/` paths**, not the upstream path, so the
links can be relocated.

## Reference genome
- `data/hg38.fa` (3.3 GB) — human reference, primary assembly + alts
- `data/hg38.fa.fai` — fasta index for fast random access via `pyfaidx`
- `data/hg38.chrom.sizes` — chrom name → length

## Regulatory elements
- `data/ENCODE_cCREs_v3.bed` — 1,063,878 candidate cis-regulatory elements
  - Columns: chrom, start, end, id1, id2, label (col 6 = e.g. "PLS",
    "pELS,CTCF-bound", "dELS", "CTCF-only,CTCF-bound", "DNase-H3K4me3")
  - Distribution: 510k dELS (no CTCF) / 278k dELS+CTCF / 172k pELS+/− CTCF /
    41k PLS+/− CTCF / 36k CTCF-only / 26k DNase-H3K4me3+/− CTCF
- `data/ccre_motif_scores.npz` — precomputed motif scores for cCREs
  (haven't explored yet)

## Promoters / TSS
- `data/gencode_v46_pc_TSS.bed` — protein-coding TSS positions
- `data/gencode_v46_all_TSS.bed` — all TSS positions

## Motifs
- `data/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt` — 2346 PFMs
  - File format: `>MA{id}\t{name}\n` then `A  [ v1 v2 ... ]\n` per base ACGT
  - Median length 9, range 4-33
  - Mixes vertebrate/plant/fungi/insect motifs (~2/3 vertebrate)
  - Use `utils.jaspar.parse_jaspar()` to load

## Helpers
- `utils/seqlib.py`: `extract()`, `extract_from_interval()`, `is_clean()`,
  `write_sequences()` — handles fasta extraction with N filtering
- `utils/jaspar.py`: `parse_jaspar()`, `sample_motif_instance()`, `consensus()`

## What's missing
- DHS Index (Meuleman et al. 2020) — would need to download from
  https://www.meuleman.org/research/dhsindex/ if needed
- ChIP-seq peaks per cell line — would need from ENCODE portal
- Per-cell-type accessibility tracks — could fetch bigwigs from ENCODE
