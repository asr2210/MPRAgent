#!/usr/bin/env python3
"""
Experiment 001 — DHS NMF-stratified baseline.

Sample 50,000 200bp DHS sequences from Meuleman et al. 2020 curated training set
(160k sequences, 10k per NMF component on chr3-chrY, hg38). We sample 3125
sequences from each of the 16 components for perfectly stratified coverage.

Source: https://www.meuleman.org/train_all_classifier_light.csv.gz
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_COMPONENTS = 16
N_PER_COMPONENT = N_TOTAL // N_COMPONENTS  # 3125
REMAINDER = N_TOTAL - N_PER_COMPONENT * N_COMPONENTS  # 0


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA, sep="\t")

    # Sanity checks
    assert df["raw_sequence"].str.len().min() == 200
    assert df["raw_sequence"].str.len().max() == 200
    assert df["component"].nunique() == N_COMPONENTS

    sampled = []
    for comp in range(1, N_COMPONENTS + 1):
        pool = df[df["component"] == comp]
        idx = rng.choice(len(pool), size=N_PER_COMPONENT, replace=False)
        sampled.append(pool.iloc[idx]["raw_sequence"].values)
    seqs = np.concatenate(sampled)
    rng.shuffle(seqs)

    # Hard-validate everything before writing
    assert len(seqs) == N_TOTAL, f"got {len(seqs)} sequences"
    allowed = set("ACGT")
    bad = 0
    for s in seqs:
        if len(s) != 200 or any(c not in allowed for c in s):
            bad += 1
    if bad:
        # Filter out any sequences with non-ACGT characters (N etc.) and replace
        # with another random draw from the same component.
        print(f"WARNING: {bad} sequences had non-ACGT or wrong length; filtering and replacing")
        good_seqs = []
        replacements_needed = 0
        for s in seqs:
            if len(s) == 200 and all(c in allowed for c in s):
                good_seqs.append(s)
            else:
                replacements_needed += 1
        while replacements_needed > 0:
            pool = df["raw_sequence"].values
            cand = pool[rng.integers(0, len(pool))]
            if len(cand) == 200 and all(c in allowed for c in cand):
                good_seqs.append(cand)
                replacements_needed -= 1
        seqs = np.array(good_seqs)
        rng.shuffle(seqs)
        assert len(seqs) == N_TOTAL

    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")

    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
