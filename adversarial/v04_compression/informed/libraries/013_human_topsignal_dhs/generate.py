#!/usr/bin/env python3
"""
013_human_topsignal_dhs — 40k human random tiles + 10k DHS sites with
highest mean_signal across 733 Meuleman biosamples. Targets
"universally active" elements as a small additive supplement.

Hypothesis: top-mean_signal DHS sites are the most universally active
regulatory elements across cell types. Mixing them into a broad random
genome at 20% adds verified cell-type-transferable signal without the
dilution disaster of 50/50 mixes (E4 = 0.493, E12 = 0.484).

Top-1% threshold (>3.70 signal) yields 35,919 sites — sample 10k uniformly.

Predicted outcomes:
- ≥0.51: high-mean_signal injection adds verified cell-type-transfer info.
- ~0.50: neutral.
- <0.49: even small DHS injection hurts; stop adding curated content.
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N_RANDOM = 40_000
N_DHS = 10_000
L = 200
HUMAN_CHROMS = ["chr1", "chr17", "chr19", "chr22"]
DATA = Path(__file__).resolve().parents[2] / "data"


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def sample_random_tiles(rng, chrom_seqs, n, L):
    lengths = np.array([len(chrom_seqs[c]) for c in chrom_seqs], dtype=float)
    weights = lengths / lengths.sum()
    chroms_list = list(chrom_seqs.keys())
    valid = set("ACGT")
    out = []
    attempts = 0
    while len(out) < n and attempts < n * 10:
        attempts += 1
        c = chroms_list[rng.choice(len(chroms_list), p=weights)]
        seq = chrom_seqs[c]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start : start + L]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    return out, attempts


def main():
    rng = np.random.default_rng(SEED)

    print("loading DHS index...")
    df = pd.read_csv(DATA / "dhs_index.txt.gz", sep="\t",
                     low_memory=False)
    std = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
    df = df[df["seqname"].isin(std)].reset_index(drop=True)
    # Top 1% by mean_signal
    threshold = df["mean_signal"].quantile(0.99)
    top = df[df["mean_signal"] >= threshold].reset_index(drop=True)
    print(f"  total: {len(df):,}, top-1% (signal≥{threshold:.2f}): {len(top):,}")

    # Sample 12k from top-1% (oversample for edge rejections)
    sample_idx = rng.choice(len(top), size=12_000, replace=False)
    sample = top.iloc[sample_idx].reset_index(drop=True)
    chroms_needed = set(sample["seqname"].unique())
    # Also need all human chroms for random tiles
    chroms_needed.update(HUMAN_CHROMS)
    print(f"loading {len(chroms_needed)} chromosomes...")
    chrom_seqs_all = {c: load_chrom(DATA / f"{c}.fa.gz") for c in chroms_needed}

    valid = set("ACGT")
    dhs_seqs = []
    for _, row in sample.iterrows():
        if len(dhs_seqs) >= N_DHS:
            break
        summit = int(row["summit"])
        start = summit - L // 2
        end = start + L
        seq = chrom_seqs_all[row["seqname"]]
        if start < 0 or end > len(seq):
            continue
        w = seq[start:end]
        if not set(w).issubset(valid):
            continue
        dhs_seqs.append(w)
    assert len(dhs_seqs) == N_DHS, len(dhs_seqs)
    print(f"DHS top-signal: {len(dhs_seqs)} sequences")

    # Random tiles from human gene-rich chroms only
    human_chrom_seqs = {c: chrom_seqs_all[c] for c in HUMAN_CHROMS}
    rand_seqs, a = sample_random_tiles(rng, human_chrom_seqs, N_RANDOM, L)
    print(f"random tiles: {len(rand_seqs)} of {a} attempts")

    out = rand_seqs + dhs_seqs
    assert len(out) == N_RANDOM + N_DHS
    rng.shuffle(out)

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")
    print(f"wrote {len(out)} → {out_path}")


if __name__ == "__main__":
    main()
