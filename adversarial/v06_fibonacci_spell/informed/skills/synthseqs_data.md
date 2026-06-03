# SynthSeqs (Meuleman 2020) — pre-curated 200bp DHS sequences

## What it is
Pre-extracted 200bp DNA sequences from the Meuleman 2020 DHS index, curated
into a balanced training set for their SynthSeqs generative model project.

- **train**: 160,000 sequences (10,000 per NMF component × 16 components)
- **val/test**: 16,000 each (1k per component × 16)
- Each sequence is exactly 200bp from hg38
- Sequences chosen for **strong topic dominance** (high `proportion` value)

## Why useful
- Avoids needing to download hg38 (3GB) + extract sequences myself
- Already balanced per NMF chromatin program → easy stratified sampling
- Includes per-sequence metadata: `total_signal` (sum across biosamples),
  `numsamples` (how many biosamples it's accessible in), `proportion`
  (purity of dominant topic)

## Download
```bash
mkdir -p data
wget -q -O data/train_synthseqs.csv.gz https://www.meuleman.org/train_all_classifier_light.csv.gz
wget -q -O data/val_synthseqs.csv.gz https://www.meuleman.org/validation_all_classifier_light.csv.gz
wget -q -O data/test_synthseqs.csv.gz https://www.meuleman.org/test_all_classifier_light.csv.gz
```

## File format
Tab-separated (despite .csv extension). Columns:
`seqname  start  end  DHS_width  summit  total_signal  numsamples  raw_sequence  component  proportion`

- `component`: integer 1–16 (NMF chromatin program; see Meuleman 2020 Fig 3)
- `proportion`: 0–1, how much of the DHS signal comes from the dominant component
- `total_signal`: sum of DHS signal across all 733 biosamples
- `numsamples`: how many biosamples this DHS is accessible in

## NMF component meanings (Meuleman 2020)
1: Placental/trophoblast, 2: Lymphoid, 3: Myeloid, 4: Erythroid,
5: Cardiac, 6: Musculoskeletal, 7: Vascular/endothelial, 8: Primitive/embryonic,
9: Neural, 10: Digestive, 11: Stromal A, 12: Renal/cancer, 13: Cancer/epithelial,
14: Pulmonary devel., 15: Organ devel./renal, 16: Tissue invariant
(Numbering may differ — confirm against the train file metadata.)

## Load in Python
```python
import pandas as pd
df = pd.read_csv("data/train_synthseqs.csv.gz", sep="\t", compression="gzip")
# df.shape = (160000, 10)
# df['raw_sequence'] is uppercase ACGT, exactly 200bp
```

## Notes
- All sequences guaranteed 200bp, no Ns observed (need to verify)
- Topic membership in synthseqs is filtered for purity — different from
  uniform DHS sampling. May behave differently than `dhs_stratified` baseline.
