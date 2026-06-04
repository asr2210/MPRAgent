#!/usr/bin/env python3
"""
Experiment 014 — GC-stratified Dirichlet(0.3).

Generate 200k Dirichlet(0.3) compositions. Compute expected GC = p[C] + p[G] for each.
Bin into 10 GC-quantile bins. Sample 5k compositions per bin → 50k. Each composition
generates one 200bp sequence.

Tests if forcing uniform GC coverage across Dirichlet(0.3) compositions improves over
random Dirichlet(0.3) draws (exp 004 = 0.0786).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_OVERSAMPLE = 200_000
N_BINS = 10
N_PER_BIN = N_TOTAL // N_BINS
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    # Step 1: oversample Dirichlet(0.3) compositions
    comps = rng.dirichlet([ALPHA] * 4, size=N_OVERSAMPLE)  # shape (N, 4) order A,C,G,T
    gc = comps[:, 1] + comps[:, 2]
    # GC quantile bins
    bins = pd.qcut(gc, q=N_BINS, labels=False, duplicates="drop")
    print(f"GC range: {gc.min():.3f}-{gc.max():.3f}, median={np.median(gc):.3f}")

    selected_comps = []
    for b in range(N_BINS):
        mask = bins == b
        candidates = np.where(mask)[0]
        if len(candidates) < N_PER_BIN:
            print(f"WARNING: bin {b} has only {len(candidates)}, sampling with replacement")
            idx = rng.choice(candidates, size=N_PER_BIN, replace=True)
        else:
            idx = rng.choice(candidates, size=N_PER_BIN, replace=False)
        selected_comps.append(comps[idx])
    selected = np.vstack(selected_comps)
    assert selected.shape == (N_TOTAL, 4)
    print(f"selected {selected.shape[0]} comps, GC range {(selected[:,1]+selected[:,2]).min():.3f}-{(selected[:,1]+selected[:,2]).max():.3f}")

    # Step 2: generate one sequence per composition
    seqs = []
    for p in selected:
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        seqs.append("".join(ALPHABET[idx]))
    seqs = np.array(seqs)
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(seqs))}")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
