#!/usr/bin/env python3
"""
Experiment 011 — Pure top-signal DHS at 50k.

Sort Meuleman 160k DHS by total_signal. Take top 50,000 (signal cutoff ≈ 11.4).
Tests whether biological signal-selection is a stronger lever than composition
variance. If yes, eval_01 should exceed the ~0.078 composition ceiling.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA, sep="\t")
    df = df.sort_values("total_signal", ascending=False).head(N_TOTAL).copy()
    print(f"selected top {len(df)} DHS by total_signal")
    print(f"signal range: {df['total_signal'].min():.2f} - {df['total_signal'].max():.2f}")
    print(f"signal median: {df['total_signal'].median():.2f}")

    seqs = df["raw_sequence"].values.copy()
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != 200 or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    n_unique = len(set(seqs))
    print(f"unique: {n_unique}")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
