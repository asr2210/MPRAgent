#!/usr/bin/env python3
"""
015_chr19_22_only — 50k random 200bp tiles from human chr19 + chr22
only. The two most gene-dense chromosomes.

chr19: ~60 Mb, ~1,400 genes (23 genes/Mb — densest in hg38)
chr22: ~51 Mb, ~740 genes  (14.5 genes/Mb)
Total: ~111 Mb. 50k tiles × 200bp = 10 Mb = ~9% coverage.
Some tile overlap expected, but not crippling.

E2 used chr1/17/19/22 (gene-dense average ~12 genes/Mb): 0.4992.
This restricts further to chr19+22 (avg ~19 genes/Mb).

Hypothesis: if gene-richness monotonically helps, this beats E2.
If narrowing hurts even when gene-rich, this is worse.

Outcomes:
- >0.51: gene-density is the right axis to push, monotonic.
- 0.49-0.51: saturated; chr1/17/19/22 already captured the gene-rich
  benefit.
- <0.49: narrowness penalty wins over gene-density gain.

seed=0. Length-weighted between chr19 and chr22.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
CHROMS = ["chr19", "chr22"]
DATA = Path(__file__).resolve().parents[2] / "data"


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}
    lengths = {c: len(s) for c, s in chrom_seqs.items()}
    print(f"chr19: {lengths['chr19']:,}; chr22: {lengths['chr22']:,}")
    weights = np.array([lengths[c] for c in CHROMS], dtype=float)
    weights /= weights.sum()
    print(f"sampling weights: {dict(zip(CHROMS, weights))}")

    valid = set("ACGT")
    out = []
    attempts = 0
    while len(out) < N and attempts < N * 10:
        attempts += 1
        c = CHROMS[rng.choice(len(CHROMS), p=weights)]
        seq = chrom_seqs[c]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start:start + L]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    print(f"wrote {len(out)} of {attempts} attempts "
          f"(reject {1 - len(out)/attempts:.3f})")
    assert len(out) == N

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
