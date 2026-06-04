#!/usr/bin/env python3
"""
Experiment 006 — 25k Dirichlet(0.3) + 25k high-signal DHS.

Tests whether high-signal biology adds to the winning composition strategy.
- Dirichlet(0.3) alone: eval_01 = 0.0786 (exp 004)
- DHS stratified alone: eval_01 = 0.0739 (exp 001)
- DHS + Dirichlet(1.0): eval_01 = 0.0765 (exp 002)

If biology adds on top of composition diversity, we should beat 0.0786.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_DIRICHLET = 25_000
N_DHS = 25_000
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def dirichlet_block(rng):
    out = []
    for _ in range(N_DIRICHLET):
        p = rng.dirichlet([ALPHA] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        out.append("".join(ALPHABET[idx]))
    return np.array(out, dtype=object)


def dhs_block(rng):
    df = pd.read_csv(DATA, sep="\t")
    # Take TOP 25% by total_signal — strong, broadly accessible DHSs.
    threshold = df["total_signal"].quantile(0.75)
    high = df[df["total_signal"] >= threshold]
    print(f"DHS pool size after top-25%-signal filter: {len(high)}")
    idx = rng.choice(len(high), size=N_DHS, replace=False)
    return high.iloc[idx]["raw_sequence"].values


def main():
    rng = np.random.default_rng(SEED)
    dir_ = dirichlet_block(rng)
    dhs = dhs_block(rng)
    seqs = np.concatenate([dir_, dhs])
    rng.shuffle(seqs)
    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
