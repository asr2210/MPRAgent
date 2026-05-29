# Skill: hg38 reference genome — download and use

## Where to download
UCSC per-chromosome FASTAs, soft-masked:
```
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr{N}.fa.gz
```
- Soft-masked = repeat regions are in lowercase, non-repeat in uppercase.
- Each chrN.fa.gz is ~10–30 MB compressed, ~50–250 MB uncompressed.
- Already downloaded in `data/`: chr19, chr20, chr21, chr22.

## Yields after filtering (drop windows with N or > 50% soft-masked, 200bp non-overlap)
| chrom | bp        | clean 200bp windows |
|-------|-----------|--------------------|
| chr19 | 58.6 M    | 114,577            |
| chr20 | 64.4 M    | 146,629            |
| chr21 | 46.7 M    | 96,684             |
| chr22 | 50.8 M    | 89,849             |
| total | 220.5 M   | **447,739**        |

So chr19–22 alone gives ~3× headroom over 150k.

For broader genome sampling, download larger chromosomes — chr1 (~250 Mb)
would alone give ~600k+ clean windows.

## Reading FASTA without pyfaidx
For single-record per-chromosome FASTAs (UCSC format), simply skip the
header line and concatenate the remaining lines. See
`libraries/002_hg38_tiled/generate.py::read_fasta_seq`.

## Default filter for "clean genomic 200bp windows"
- exclude any window containing chars outside `ACGTacgt` (drops N tracts at
  telomeres / centromeres / gaps)
- exclude windows with > 50% lowercase chars (drops the most heavily
  repeat-masked regions; keeps some repeats so we don't lose all SINE/LINE)
- upper-case the survivors before writing
- sample without replacement with a fixed seed
