#!/usr/bin/env python3
"""
Experiment 019 — 4D-stratified Dirichlet(0.3).

Generate 500k Dirichlet(0.3) compositions. Bin each base's proportion into 5
quantile bins → 5×5×5×5 = 625 4D bins. Sample 80 compositions per bin (with
replacement if needed) → 50,000.

Tests if forcing uniform 4D composition coverage breaks the Dirichlet(0.3) ceiling.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_OVERSAMPLE = 500_000
N_TOTAL = 50_000
N_BINS_PER_AXIS = 5
N_PER_BIN = 80  # 80 * 625 = 50000
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    comps = rng.dirichlet([ALPHA] * 4, size=N_OVERSAMPLE)  # shape (N, 4) order A,C,G,T

    # Define bin edges from quantiles of each marginal
    edges = []
    for k in range(4):
        e = np.quantile(comps[:, k], np.linspace(0, 1, N_BINS_PER_AXIS + 1))
        e[0] = -np.inf
        e[-1] = np.inf
        edges.append(e)

    bins = np.zeros((N_OVERSAMPLE, 4), dtype=int)
    for k in range(4):
        bins[:, k] = np.searchsorted(edges[k][1:-1], comps[:, k])
    # bin id = b0*125 + b1*25 + b2*5 + b3
    bin_id = bins[:, 0] * 125 + bins[:, 1] * 25 + bins[:, 2] * 5 + bins[:, 3]

    sel_idx = []
    n_bins = N_BINS_PER_AXIS ** 4
    empty_bins = 0
    for b in range(n_bins):
        cand = np.where(bin_id == b)[0]
        if len(cand) == 0:
            # empty bin: skip (will sample extra from non-empty)
            empty_bins += 1
            continue
        if len(cand) >= N_PER_BIN:
            sel_idx.extend(rng.choice(cand, size=N_PER_BIN, replace=False))
        else:
            sel_idx.extend(cand)
            short = N_PER_BIN - len(cand)
            sel_idx.extend(rng.choice(cand, size=short, replace=True))
    print(f"empty bins: {empty_bins}/{n_bins}")

    # Top up to 50000 if we lost any due to empty bins (no empty since searchsorted always assigns)
    while len(sel_idx) < N_TOTAL:
        sel_idx.append(int(rng.integers(N_OVERSAMPLE)))
    sel_idx = sel_idx[:N_TOTAL]

    selected = comps[sel_idx]
    print(f"selected {selected.shape[0]} comps")

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
