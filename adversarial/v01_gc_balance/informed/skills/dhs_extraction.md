# Skill: DHS region extraction (Meuleman 2020 → hg38 200bp sequences)

**Data**:
- `data/DHS_Index_hg38.txt.gz` (~3.59M DHSs, ~90MB gzipped)
- `data/hg38.fa` (3.2GB decompressed)
- `data/DHS_Index_and_Vocabulary_metadata.tsv` (733 biosamples)

**Index format** (TSV with header):
| col | name | meaning |
|---|---|---|
| 1 | seqname | chromosome (chr1..chr22, chrX, chrY) |
| 2 | start | DHS start (hg38) |
| 3 | end | DHS end |
| 4 | identifier | string id |
| 5 | mean_signal | average DNase signal across biosamples in which DHS is called |
| 6 | numsamples | number of biosamples (of 733) in which this DHS is called |
| 7 | summit | summit position (peak of accessibility) |
| 8 | core_start | core motif start |
| 9 | core_end | core motif end |
| 10 | component | dominant NMF topic (16 named topics) |

**NMF topics and counts** (total 3.59M):
```
Primitive / embryonic     626541
Neural                    461478
Stromal B                 404883
Lymphoid                  280192
Placental / trophoblast   264980
Musculoskeletal           216894
Cancer / epithelial       188489
Myeloid / erythroid       186616
Organ devel. / renal      159137
Tissue invariant          157670
Digestive                 144681
Renal / cancer            144087
Cardiac                   118869
Pulmonary devel.           96369
Vascular / endothelial     84826
Stromal A                  56186
```
**K562 ≈ Myeloid/erythroid**, **HepG2 ≈ Digestive (hepatic)**, **SK-N-SH ≈ Neural**.

**Pattern for sampling 50k 200bp sequences centred on DHS summit**:
```python
import gzip
import pandas as pd
import numpy as np
from pyfaidx import Fasta

df = pd.read_csv('data/DHS_Index_hg38.txt.gz', sep='\t', compression='gzip')
fa = Fasta('data/hg38.fa', as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df['seqname'].unique()}

# Filter to DHSs where summit ± 100 fits within chromosome
half = 100
df['win_start'] = df['summit'] - half
df['win_end'] = df['summit'] + half
df = df[(df['win_start'] >= 0) & (df.apply(lambda r: r['win_end'] <= chrom_len[r['seqname']], axis=1))]

# Sample 50k (with or without replacement)
sample = df.sample(n=50_000, replace=False, random_state=42)

def get_seq(chrom, s, e):
    return str(fa[chrom][s:e]).upper()

seqs = [get_seq(r.seqname, r.win_start, r.win_end) for r in sample.itertuples()]
# Drop any with non-ACGT (e.g. N) — re-sample replacements
```

**Edge cases**:
- `pyfaidx` returns lowercase or uppercase depending on `sequence_always_upper`. Force upper.
- Genome contains `N` bases (gaps, centromeres, telomeres). DHSs are usually in mapped regions but a few may have Ns. Filter or replace.
- Some DHSs are very small (60bp); using summit ± 100 always yields exactly 200bp.
- pyfaidx is fast for random access after `.fai` index built (one-time ~30s).

**Sampling strategies** (built up over experiments):
- `uniform` — sample any DHS with equal prob. Baseline ≈ `dhs_random` (eval_01 ≈ 0.71).
- `topic_proportional` — sample proportional to NMF topic size (= uniform, since each region has one max-topic label).
- `topic_balanced` — equal counts per topic (50k / 16 ≈ 3125 per topic). Baseline ≈ `dhs_stratified` (eval_01 ≈ 0.71).
- `numsamples_weighted` — sample with weight ∝ numsamples (broadly-active DHSs). NEW, untested.
- `mean_signal_weighted` — sample with weight ∝ mean_signal. NEW, untested.
- `cell_type_excluded` — drop regions whose top topic matches a labeling cell (K562 ~ Myeloid, HepG2 ~ Digestive, SK-N-SH ~ Neural). Tests generalisation.
