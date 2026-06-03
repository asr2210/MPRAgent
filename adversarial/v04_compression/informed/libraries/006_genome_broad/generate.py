#!/usr/bin/env python3
"""
006_genome_broad — random tiling from all 24 hg38 chromosomes.

Robustness check on 002. 002 used only chr1/17/19/22 (gene-rich
biased subset). Here we sample from chr1–chr22 + chrX + chrY,
weighted by chromosome length.

If 006 ≈ 002 (within 0.03 on eval_01): real-DNA effect is robust to
chromosome choice; lock in this design as my baseline.
If 006 >> 002: broader sampling helps further — pursue.
If 006 << 002: chr1/17/19/22 had a special advantage (gene density).
Need to identify and replicate the advantage.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        first = f.readline()
        assert first.startswith(">")
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    chrom_seqs = {}
    for c in CHROMS:
        path = DATA_DIR / f"{c}.fa.gz"
        chrom_seqs[c] = load_chrom(path)
    total = sum(len(s) for s in chrom_seqs.values())
    print(f"loaded {len(chrom_seqs)} chromosomes, total {total:,} bp")

    weights = np.array([len(chrom_seqs[c]) - L for c in CHROMS], dtype=np.float64)
    weights /= weights.sum()

    valid = set("ACGT")
    out = []
    attempts = 0
    while len(out) < N:
        batch = max(N - len(out), 1) * 2
        cidx = rng.choice(len(CHROMS), size=batch, p=weights)
        for ci in cidx:
            if len(out) >= N:
                break
            attempts += 1
            seq = chrom_seqs[CHROMS[ci]]
            start = rng.integers(0, len(seq) - L)
            w = seq[start : start + L]
            if not set(w).issubset(valid):
                continue
            out.append(w)
    print(f"wrote {len(out)} of {attempts} attempts ({1 - len(out)/attempts:.3f} reject)")

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
