# Skill: extracting 200bp sequences from the Meuleman DHS Index

## Data sources
- DHS index file:
  `https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz`
  (~87 MB gz; ~3.59M elements; header row + tab-separated)
- hg38 reference: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
  (~939 MB gz, ~3.1 GB unpacked)

Stored locally at `data/dhs_index.txt.gz` and `data/hg38.fa` (+ `.fai`).

## DHS index columns
```
seqname start end identifier mean_signal numsamples summit core_start core_end component
```

- `summit`: best single-bp summit position (0-based). Use as center for 200bp windows.
- `numsamples`: # biosamples this element is open in. Low = cell-type-specific.
- `mean_signal`: aggregated DNase signal. Higher = stronger accessibility.
- `component`: one of 16 NMF dominant components (text label, see below). The full
  16-dim NMF loadings are NOT in this file — only the dominant component label.

## 16 NMF components (counts in the index)
Primitive / embryonic 626k, Neural 461k, Stromal B 405k, Lymphoid 280k,
Placental / trophoblast 265k, Musculoskeletal 217k, Cancer / epithelial 188k,
Myeloid / erythroid 187k, Organ devel. / renal 159k, Tissue invariant 158k,
Digestive 145k, Renal / cancer 144k, Cardiac 119k, Pulmonary devel. 96k,
Vascular / endothelial 85k, Stromal A 56k.

## Recipe: extract 200bp windows centered on summit
```python
import gzip
import numpy as np
from pyfaidx import Fasta

fa = Fasta("data/hg38.fa", as_raw=False, sequence_always_upper=True)
chroms = set(fa.keys())  # filter to known canonical chroms

records = []  # (chrom, summit, mean_signal, numsamples, component)
with gzip.open("data/dhs_index.txt.gz", "rt") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, _, _, _, sig, nsamp, summit, _, _, comp = parts
        if chrom not in chroms:
            continue
        records.append((chrom, int(summit), float(sig), int(nsamp), comp))

rng = np.random.default_rng(0)
idx = rng.choice(len(records), size=N_TARGET, replace=False)
SEQ_LEN = 200
HALF = SEQ_LEN // 2
seqs = []
for i in idx:
    chrom, summit, _, _, _ = records[i]
    L = len(fa[chrom])
    start = max(0, summit - HALF)
    end = start + SEQ_LEN
    if end > L:
        end = L
        start = end - SEQ_LEN
    s = str(fa[chrom][start:end]).upper()
    if len(s) != SEQ_LEN or set(s) - set("ACGT"):
        continue  # skip N-containing or truncated windows
    seqs.append(s)
```

## Sampling variants
- **dhs_random**: uniform `rng.choice(len(records), size=N)`
- **dhs_stratified**: group by `component`, sample N/16 per component
- **dhs_specific** (proxy for topic weighting): weight ∝ 1/numsamples (favor specific)
- **dhs_strong**: weight ∝ mean_signal (favor strong peaks)

## Gotchas
- Some DHS sites overlap N stretches in hg38 → check that extracted sequence
  contains only ACGT.
- ~10% of the index is on alt/random contigs — drop these by checking
  `chrom in fa.keys()` and filtering to canonical chr1..22,X,Y,M.
- Use `pyfaidx.Fasta(... as_raw=False, sequence_always_upper=True)` for case
  normalization; or `.seq.upper()` post hoc.
