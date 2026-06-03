#!/usr/bin/env python3
"""
026_dhs_density_weighted — same scheme as E24 (250kb bins, squared
weighting), but using DHS site density (Meuleman index) instead of
gene density as the per-bin weight.

DHS density may match the eval distribution better than gene density
because DHS sites are by definition "accessible/regulatory" and the
eval is regulatory (MPRA activity).

Compare:
- E24 (gene-density squared, 250kb): 0.5084
- E26 (DHS-density squared, 250kb): ?

If E26 > E24: DHS-density is a better proxy. Push further.
If E26 ≈ E24: similar signal; gene-rich and DHS-dense regions
overlap.
If E26 < E24: gene-density is the dominant useful signal.
"""
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

N = 50_000
L = 200
BIN_SIZE = 250_000
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
    print("loading DHS index...")
    dhs = pd.read_csv(DATA / "dhs_index.txt.gz", sep="\t", low_memory=False)
    dhs = dhs[dhs["seqname"].isin(set(CHROMS))].reset_index(drop=True)
    print(f"  {len(dhs):,} DHS sites on standard chroms")

    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}

    bin_weights = []
    for c in CHROMS:
        chrom_len = len(chrom_seqs[c])
        n_bins = chrom_len // BIN_SIZE
        counts = np.zeros(n_bins, dtype=np.int32)
        sites = dhs[dhs["seqname"] == c]["summit"].values
        for s in sites:
            b = int(s) // BIN_SIZE
            if 0 <= b < n_bins:
                counts[b] += 1
        for b in range(n_bins):
            bin_weights.append((c, b, (counts[b] + EPS) ** EXP))

    weights = np.array([w for _, _, w in bin_weights], dtype=float)
    weights /= weights.sum()
    print(f"{len(bin_weights):,} bins; top-50 mass: {sum(np.sort(weights)[-50:]):.4f}")

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
