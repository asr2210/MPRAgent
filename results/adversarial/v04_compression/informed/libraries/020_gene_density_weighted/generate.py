#!/usr/bin/env python3
"""
020_gene_density_weighted — 50k random 200bp tiles sampled across
all 24 hg38 chromosomes, but weighted by LOCAL gene density per 1Mb bin.

Computes gene count per 1Mb bin from RefSeq GTF, samples bin
proportional to (gene_count + epsilon), then samples uniformly within
bin. Tests if continuous gene-density weighting beats hard chromosome
filter (E2 chr1/17/19/22 = 0.4992; E6 all-chroms = 0.4823).

Expected: between 0.49 and 0.50. Probably saturated. If >0.50, gene
density is a continuously-monotonic axis; if <0.49, the weighting
schema introduces some unwanted bias.

seed=0.
"""
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
BIN_SIZE = 1_000_000  # 1 Mb
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
DATA = Path(__file__).resolve().parents[2] / "data"
EPS = 0.5  # smooth: tiles per bin > 0 even for gene-empty bins


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    # Parse RefSeq gene starts
    print("parsing GTF for gene starts...")
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
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            t = start if strand == "+" else end
            key = (chrom, t)
            if key in seen:
                continue
            seen.add(key)
            gene_starts[chrom].append(t)
    total_genes = sum(len(v) for v in gene_starts.values())
    print(f"  {total_genes:,} unique transcript starts")

    print("loading chromosomes...")
    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}

    # Build per-bin gene counts
    bin_weights = []  # list of (chrom, bin_index, weight)
    for c in CHROMS:
        chrom_len = len(chrom_seqs[c])
        n_bins = chrom_len // BIN_SIZE
        counts = np.zeros(n_bins, dtype=np.int32)
        for g in gene_starts.get(c, []):
            b = g // BIN_SIZE
            if 0 <= b < n_bins:
                counts[b] += 1
        for b in range(n_bins):
            bin_weights.append((c, b, counts[b] + EPS))

    weights = np.array([w for _, _, w in bin_weights], dtype=float)
    weights /= weights.sum()
    print(f"  {len(bin_weights):,} bins; top10 most-dense:")
    top_idx = np.argsort(weights)[-10:][::-1]
    for i in top_idx:
        c, b, w = bin_weights[i]
        print(f"    {c}:{b*BIN_SIZE//1_000_000}Mb weight={w:.1f}")

    valid = set("ACGT")
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

    print(f"wrote {len(out)} of {attempts} attempts "
          f"(reject {1 - len(out)/attempts:.3f})")
    assert len(out) == N

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
