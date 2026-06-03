#!/usr/bin/env python3
"""
027_combined_density — combine gene + cCRE density signals.

Per 250kb bin, compute:
  combined = z(gene_count) + z(ccre_count)
  weight = max(combined, 0) + EPS, then squared.

z(x) = (x - mean) / std for non-zero bins.

E24 (gene^2, 250kb):  0.5084
E26 (DHS^2, 250kb):   0.5067
E27 (combined):       ?

If E27 > E24: combining signals beats either alone.
If E27 ≈ E24: gene-density is dominant.
If E27 < E24: combining introduces noise / over-narrowing.

3-seed.
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
    # Gene starts
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
    print(f"transcripts: {sum(len(v) for v in gene_starts.values()):,}")

    # cCRE midpoints
    print("loading cCRE...")
    ccre = pd.read_csv(DATA / "cCRE_hg38.bed", sep="\t", header=None,
                       names=["chrom", "start", "end", "dccid", "sccid", "class"])
    ccre["mid"] = (ccre["start"] + ccre["end"]) // 2
    ccre = ccre[ccre["chrom"].isin(set(CHROMS))].reset_index(drop=True)
    print(f"  cCREs: {len(ccre):,}")

    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}

    # Bin counts
    bin_meta = []  # (chrom, b, gene_count, ccre_count)
    for c in CHROMS:
        chrom_len = len(chrom_seqs[c])
        n_bins = chrom_len // BIN_SIZE
        gcounts = np.zeros(n_bins, dtype=np.int32)
        ccounts = np.zeros(n_bins, dtype=np.int32)
        for g in gene_starts.get(c, []):
            b = g // BIN_SIZE
            if 0 <= b < n_bins:
                gcounts[b] += 1
        sub = ccre[ccre["chrom"] == c]["mid"].values
        for s in sub:
            b = int(s) // BIN_SIZE
            if 0 <= b < n_bins:
                ccounts[b] += 1
        for b in range(n_bins):
            bin_meta.append((c, b, gcounts[b], ccounts[b]))

    gcs = np.array([m[2] for m in bin_meta], dtype=float)
    ccs = np.array([m[3] for m in bin_meta], dtype=float)
    gz = (gcs - gcs.mean()) / (gcs.std() + 1e-9)
    cz = (ccs - ccs.mean()) / (ccs.std() + 1e-9)
    combined = np.maximum(gz + cz, 0) + EPS
    weights = combined ** EXP
    weights /= weights.sum()
    print(f"{len(bin_meta):,} bins; top-50 mass: {sum(np.sort(weights)[-50:]):.4f}")

    valid = set("ACGT")
    out_dir = Path(__file__).resolve().parent
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        out = []
        attempts = 0
        while len(out) < N and attempts < N * 10:
            attempts += 1
            bi = rng.choice(len(bin_meta), p=weights)
            c, b, _, _ = bin_meta[bi]
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
