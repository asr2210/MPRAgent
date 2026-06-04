#!/usr/bin/env python3
"""
Experiment 016 — 25k Dirichlet(0.3) + 25k shuffled DHS.

Combine best synthetic (Dirichlet(0.3), exp 004 = 0.0786) with best biology
(shuffled DHS, exp 013 = 0.0754).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_HALF = N_TOTAL // 2
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)

    # Half 1: Dirichlet(0.3) synthetic
    dirichlet_seqs = []
    for _ in range(N_HALF):
        p = rng.dirichlet([ALPHA] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        dirichlet_seqs.append("".join(ALPHABET[idx]))

    # Half 2: shuffled DHS NMF-stratified (16 components × ~1562)
    df = pd.read_csv(DATA, sep="\t")
    components = sorted(df["component"].unique())
    n_per_component = N_HALF // len(components)  # ~1562
    extra = N_HALF - n_per_component * len(components)
    bio_seqs = []
    for i, c in enumerate(components):
        n_take = n_per_component + (1 if i < extra else 0)
        pool = df[df["component"] == c]
        idx = rng.choice(len(pool), size=n_take, replace=False)
        bio_seqs.extend(pool.iloc[idx]["raw_sequence"].values)
    assert len(bio_seqs) == N_HALF

    # Shuffle each bio sequence internally
    shuffled_bio = []
    for s in bio_seqs:
        arr = np.array(list(s))
        rng.shuffle(arr)
        shuffled_bio.append("".join(arr))

    all_seqs = np.array(dirichlet_seqs + shuffled_bio)
    rng.shuffle(all_seqs)

    assert len(all_seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in all_seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(all_seqs))}")
    with open(OUT, "w") as fh:
        for s in all_seqs:
            fh.write(s + "\n")
    print(f"wrote {len(all_seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
