#!/usr/bin/env python3
"""
Experiment 008 — DHS sequences GC-content-stratified.

Take Meuleman 160k DHS pool. Compute GC content per sequence. Bin into 10
quantile bins of equal pool size. Sample 5,000 per bin → 50,000 sequences.

Tests whether biology (real DHS) becomes useful when EXPLICITLY arranged with
compositional spread, vs random/NMF-stratified DHS where compositions are
concentrated near 50% GC.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_BINS = 10
N_PER_BIN = N_TOTAL // N_BINS  # 5000


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA, sep="\t")
    # Compute GC content per sequence.
    df["gc"] = df["raw_sequence"].apply(lambda s: (s.count("G") + s.count("C")) / len(s))
    # Quantile bins
    df["gc_bin"] = pd.qcut(df["gc"], q=N_BINS, labels=False, duplicates="drop")
    print(f"DHS GC content: min={df['gc'].min():.3f}, max={df['gc'].max():.3f}, mean={df['gc'].mean():.3f}, std={df['gc'].std():.3f}")
    print("bin counts:", df["gc_bin"].value_counts().sort_index().tolist())

    seqs = []
    for b in range(N_BINS):
        pool = df[df["gc_bin"] == b]
        idx = rng.choice(len(pool), size=N_PER_BIN, replace=len(pool) < N_PER_BIN)
        seqs.extend(pool.iloc[idx]["raw_sequence"].values)
    seqs = np.array(seqs)
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != 200 or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
