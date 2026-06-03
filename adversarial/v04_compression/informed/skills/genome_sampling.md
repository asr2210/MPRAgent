# Skill: Sampling 200bp windows from hg38

## When to use
Any time you need to draw random or coordinate-targeted 200bp windows
from human reference DNA.

## Data location
- Per-chromosome gzipped fasta: `data/chr{1..22,X,Y}.fa.gz`
- Each file has 1 header line then sequence in 50-char chunks.
- All standard nucleotides are uppercase A/C/G/T, with masked / repeat
  regions appearing as lowercase or N. The Markov experiment showed
  that loading with `.upper()` blends repeats into the run; if you want
  to keep repeat awareness, don't upper-case (and filter on `set("ACGT")`
  *before* upper-casing).

## Loader pattern
```python
import gzip
from pathlib import Path

def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        first = f.readline()
        assert first.startswith(">")
        return "".join(line.strip() for line in f).upper()
```

## Uniform random tile sampler
- Reject windows containing any non-ACGT base.
- Weight per-chromosome sampling by `len(chrom) - L`.
- ~7% reject rate on chr1/17/19/22 (mostly N runs in centromeres/telomeres).
- Sampling 50k tiles + reject-loop takes ~10s end-to-end on chr1+17+19+22.

```python
def sample_tiles(rng, chrom_seqs, n, L=200):
    weights = np.array([len(s) - L for s in chrom_seqs.values()], dtype=float)
    weights /= weights.sum()
    chroms = list(chrom_seqs)
    valid = set("ACGT")
    out = []
    while len(out) < n:
        ci = rng.choice(len(chroms), size=2 * (n - len(out)), p=weights)
        for i in ci:
            if len(out) >= n: break
            seq = chrom_seqs[chroms[i]]
            start = rng.integers(0, len(seq) - L)
            w = seq[start:start+L]
            if set(w).issubset(valid):
                out.append(w)
    return out
```

## Coordinate-targeted extraction
For BED-style targets (e.g. DHS index):
```python
df = pd.read_csv('data/dhs_index.txt.gz', sep='\t')
# columns: seqname, start, end, identifier, mean_signal, numsamples,
#          summit, core_start, core_end, component
# extract 200bp centered on summit:
def extract(row, chrom_seqs, L=200):
    s = row['summit'] - L // 2
    return chrom_seqs[row['seqname']][s : s + L]
```

## Common pitfalls
- chr1/17/19/22 are gene-rich (high gene density). Sampling only
  from them gives a non-representative slice. For unbiased real-DNA
  experiments use all 24 chromosomes weighted by length.
- N-runs at centromeres/telomeres reject ~7% of windows; over-sample.
- chrY is tiny (~57 Mb after N-filter) and male-specific, contributes
  marginally with length-weighting.
- `total_signal` in SynthSeqs is the sum-of-DNase-signal across 733
  biosamples; `numsamples` is the count of biosamples where the DHS
  is detected — different concepts.
