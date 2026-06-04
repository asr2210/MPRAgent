#!/usr/bin/env python3
"""
014_seed_variance — same design as E2 (human chr1/17/19/22 random
200bp tiles), but generated with seeds 0, 1, 2 to test if prepare.py
runs multi-seed when multiple sequence files are present, and to
quantify single-seed variance.

Critical calibration: many of my single-seed Δ measurements are in
the ±0.005-0.01 range. If seed variance is that large, half of my
"insights" are noise.

E2 (seed 0): eval_01 = 0.4992. Repeats here with seeds 0, 1, 2.
- If n_seeds=1 still in result.json: multi-seed not auto-picked-up.
  Falls back to single seed=0 (will match E2 exactly).
- If multi-seed runs: standard error tells me my decision threshold.
"""
import gzip
from pathlib import Path

import numpy as np

N = 50_000
L = 200
CHROMS = ["chr1", "chr17", "chr19", "chr22"]
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
SEEDS = [0, 1, 2]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def sample_tiles(rng, chrom_seqs, n, L):
    chroms_list = list(chrom_seqs.keys())
    weights = np.array([len(chrom_seqs[c]) - L for c in chroms_list],
                       dtype=float)
    weights /= weights.sum()
    valid = set("ACGT")
    out = []
    attempts = 0
    while len(out) < n and attempts < n * 5:
        attempts += 1
        ci = rng.choice(len(chroms_list), p=weights)
        seq = chrom_seqs[chroms_list[ci]]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start:start + L]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    return out


def main():
    print("loading chromosomes...")
    chrom_seqs = {c: load_chrom(DATA_DIR / f"{c}.fa.gz") for c in CHROMS}
    print(f"  total {sum(len(s) for s in chrom_seqs.values()):,} bp")

    out_dir = Path(__file__).resolve().parent
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        seqs = sample_tiles(rng, chrom_seqs, N, L)
        assert len(seqs) == N
        out_path = out_dir / f"sequences_{seed}.txt"
        with open(out_path, "w") as f:
            f.write("\n".join(seqs))
            f.write("\n")
        print(f"  seed {seed}: wrote {len(seqs)} → {out_path.name}")


if __name__ == "__main__":
    main()
