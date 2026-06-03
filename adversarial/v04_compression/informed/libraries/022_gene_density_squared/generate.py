#!/usr/bin/env python3
"""
022_gene_density_squared — stronger gene-density weighting than E20/E21.

E20/E21 used weight = (gene_count + 0.5). This uses weight =
(gene_count + 0.5)^2, which more aggressively up-weights the
gene-richest bins.

Predicted:
- Monotonic up to 0.504-0.506: gene-density axis has more room.
- Saturated ~0.502: diminishing returns; E20/E21 close to optimum.
- Narrow penalty kicks in: <0.500 (means quadratic is too aggressive).

3-seed for reliable comparison vs E21 (0.5023).
"""
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np

N = 50_000
L = 200
BIN_SIZE = 1_000_000
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
DATA = Path(__file__).resolve().parents[2] / "data"
EPS = 0.5
EXP = 2.0
SEEDS = [0, 1, 2]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    seen = set()
    gene_starts = defaultdict(list)
    with gzip.open(DATA / "hg38.refGene.gtf.gz", "rt") as f:
        for line in f:
            parts = line.split("\t")
            if parts[2] != "transcript":
                continue
            chrom = parts[0]
            if chrom not in set(CHROMS):
                continue
            t = int(parts[3]) if parts[6] == "+" else int(parts[4])
            key = (chrom, t)
            if key in seen:
                continue
            seen.add(key)
            gene_starts[chrom].append(t)

    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}

    bin_weights = []
    for c in CHROMS:
        chrom_len = len(chrom_seqs[c])
        n_bins = chrom_len // BIN_SIZE
        counts = np.zeros(n_bins, dtype=np.int32)
        for g in gene_starts.get(c, []):
            b = g // BIN_SIZE
            if 0 <= b < n_bins:
                counts[b] += 1
        for b in range(n_bins):
            bin_weights.append((c, b, (counts[b] + EPS) ** EXP))

    weights = np.array([w for _, _, w in bin_weights], dtype=float)
    weights /= weights.sum()
    print(f"top-5 most-weighted bins:")
    top_idx = np.argsort(weights)[-5:][::-1]
    for i in top_idx:
        c, b, w = bin_weights[i]
        print(f"  {c}:{b}Mb (weight {w:.0f}, fraction {weights[i]:.4f})")

    valid = set("ACGT")
    out_dir = Path(__file__).resolve().parent
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        out = []
        attempts = 0
        while len(out) < N and attempts < N * 10:
            attempts += 1
            bi = rng.choice(len(bin_weights), p=weights)
            c, b, _ = bin_weights[bi]
            seq = chrom_seqs[c]
            bin_start = b * BIN_SIZE
            bin_end = min(bin_start + BIN_SIZE, len(seq))
            if bin_end - bin_start < L:
                continue
            start = int(rng.integers(bin_start, bin_end - L))
            w = seq[start:start + L]
            if not set(w).issubset(valid):
                continue
            out.append(w)
        assert len(out) == N
        out_path = out_dir / f"sequences_{seed}.txt"
        with open(out_path, "w") as f:
            f.write("\n".join(out))
            f.write("\n")
        print(f"seed {seed}: wrote {len(out)}")


if __name__ == "__main__":
    main()
