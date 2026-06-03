#!/usr/bin/env python3
"""
Experiment 013 — Internally-shuffled DHS.

Take Meuleman 160k DHS NMF-stratified (same as exp 001). Permute each sequence's
base order internally (preserves per-sequence composition exactly, destroys all
motifs and positional information).

Direct test of theory v4: if model learns only overall composition, shuffled DHS
should give the same eval as unshuffled DHS (exp 001 = 0.0739).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_PER_COMPONENT = 3125  # matches exp 001 (16 components × 3125 = 50k)
SEQ_LEN = 200


def shuffle_seq(seq, rng):
    arr = np.array(list(seq))
    rng.shuffle(arr)
    return "".join(arr)


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA, sep="\t")
    components = sorted(df["component"].unique())
    print(f"components: {components}")

    seqs = []
    for c in components:
        pool = df[df["component"] == c]
        idx = rng.choice(len(pool), size=N_PER_COMPONENT, replace=False)
        seqs.extend(pool.iloc[idx]["raw_sequence"].values)

    # Shuffle each sequence internally
    shuffled = [shuffle_seq(s, rng) for s in seqs]
    arr = np.array(shuffled)
    rng.shuffle(arr)  # also shuffle order

    assert len(arr) == 50_000
    allowed = set("ACGT")
    bad = sum(1 for s in arr if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    # verify composition preserved by spot-checking first 100
    for i in range(100):
        orig_comp = sorted(seqs[i])
        # i was reordered; can't compare element-by-element after final shuffle. skip.
    print(f"unique: {len(set(arr))}")
    with open(OUT, "w") as fh:
        for s in arr:
            fh.write(s + "\n")
    print(f"wrote {len(arr)} sequences to {OUT}")


if __name__ == "__main__":
    main()
