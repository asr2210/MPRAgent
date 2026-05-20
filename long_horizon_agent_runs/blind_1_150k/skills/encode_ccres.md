# Skill: ENCODE V3 cCREs for MPRA library design

## Source
GRCh38 candidate cis-regulatory elements (V3) from ENCODE SCREEN:
```
https://downloads.wenglab.org/V3/GRCh38-cCREs.bed
```
~52 MB BED. 1,063,878 elements. Stored at `data/encode_ccres_hg38.bed`.

## Format
6-column BED:
```
chrom  start  end  acc_d  acc_e  class
chr1   104896 105048 EH38D4327509  EH38E2776520  CTCF-only,CTCF-bound
```
`class` column has classes (comma-separated):
- `PLS` — promoter-like signature (DNase + H3K4me3, close to TSS)
- `pELS` — proximal enhancer-like (DNase + H3K27ac, <2kb from TSS)
- `dELS` — distal enhancer-like (DNase + H3K27ac, >2kb from TSS)
- `DNase-H3K4me3` — DNase + H3K4me3, no TSS overlap
- `CTCF-only` — DNase + CTCF only
- `CTCF-bound` — additional flag on the others
Some entries have multiple classes (e.g., `dELS,CTCF-bound`).

## Class distribution (hg38 V3, all cell types)
| class                    | count   | fraction |
|--------------------------|---------|----------|
| dELS                     | 510,920 | 48.0%   |
| dELS,CTCF-bound          | 278,280 | 26.2%   |
| pELS,CTCF-bound          |  96,781 |  9.1%   |
| pELS                     |  75,246 |  7.1%   |
| CTCF-only,CTCF-bound     |  35,839 |  3.4%   |
| PLS,CTCF-bound           |  31,447 |  3.0%   |
| DNase-H3K4me3            |  17,627 |  1.7%   |
| PLS                      |   9,444 |  0.9%   |
| DNase-H3K4me3,CTCF-bound |   8,294 |  0.8%   |

Distal enhancer-like (with or without CTCF) is the dominant class (74%).
Promoter-like (PLS, pELS variants) is ~20%.

## Standard window extraction (200bp centered)
```python
mid = (start + end) // 2
window = chrom_seq[mid - 100 : mid + 100]
```
After filtering (drop OOB, any non-ACGT char, >50% softmasked):
**~745,400 windows survive** out of 1,063,878 (~70%).

## Empirical performance (mean_r averaged across the 14 anonymous evals)
Reference values from this project's experiment series:
- Uniform random 200bp (no enrichment):        **0.820**
- hg38 chr19-22 random tiles:                  **0.873**
- ENCODE cCRE-centered (random sample of all): **0.886**
- 50/50 cCRE + uniform random hybrid:          **0.880**
- 50/50 cCRE + JASPAR motif-embedded hybrid:   **0.883**

So replacing uniform-random with cCREs adds ~+0.066 mean_r. Adding a 50%
non-cCRE diversity component costs ~−0.006 to −0.003 on mean (depending
on diversity source), but consistently lifts eval_08 by ~+0.024
(super-additive: 0.916 > both 0.892 cCRE and 0.908 random).

## Practical pitfalls
- Many cCREs lie in or near repeat regions (~30% of cCRE-centered windows
  fail the >50%-softmasked filter). Be aware that the kept set is biased
  away from repeat-rich regulatory regions (e.g., transposon-derived
  enhancers). For a balanced library, consider relaxing the repeat filter.
- The natural cCRE class distribution is dominated by distal enhancers
  (74%). For class diversity, sample with class quotas instead of
  uniformly.
- cCREs span multiple cell types but aren't cell-type specific in the
  bed file. SCREEN provides per-cell-type signal data if you want to
  filter to "K562-active" or "HepG2-active" subsets.
