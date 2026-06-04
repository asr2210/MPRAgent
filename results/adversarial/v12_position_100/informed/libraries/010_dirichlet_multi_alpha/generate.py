#!/usr/bin/env python3
"""
Experiment 010 — Multi-alpha Dirichlet mix.

10,000 sequences each at alpha ∈ {0.1, 0.3, 0.7, 1.5, 3.0} = 50,000 total.
Tests whether mixing variance regimes gives the model more useful coverage
of the composition→activity curve than a single alpha (exp 004 = 0.0786).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_PER_ALPHA = 10_000
ALPHAS = [0.1, 0.3, 0.7, 1.5, 3.0]
SEQ_LEN = 200
ALPHABET = np.array(list("ACGT"))

assert N_PER_ALPHA * len(ALPHAS) == 50_000


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for alpha in ALPHAS:
        for _ in range(N_PER_ALPHA):
            p = rng.dirichlet([alpha] * 4)
            idx = rng.choice(4, size=SEQ_LEN, p=p)
            seqs.append("".join(ALPHABET[idx]))
    seqs = np.array(seqs)
    rng.shuffle(seqs)

    assert len(seqs) == 50_000
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    n_unique = len(set(seqs))
    print(f"unique sequences: {n_unique} (expect ~50000)")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
