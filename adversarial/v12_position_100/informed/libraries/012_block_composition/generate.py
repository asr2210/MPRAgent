#!/usr/bin/env python3
"""
Experiment 012 — Block-composition Dirichlet(0.3).

Each 200bp sequence = two 100bp blocks. Each block has its own independent
Dirichlet(0.3) composition. Tests whether positional composition variation
adds learnable signal beyond uniform per-sequence composition.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
BLOCK_LEN = 100
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for _ in range(N_TOTAL):
        p1 = rng.dirichlet([ALPHA] * 4)
        p2 = rng.dirichlet([ALPHA] * 4)
        idx1 = rng.choice(4, size=BLOCK_LEN, p=p1)
        idx2 = rng.choice(4, size=BLOCK_LEN, p=p2)
        seq = "".join(ALPHABET[idx1]) + "".join(ALPHABET[idx2])
        seqs.append(seq)
    seqs = np.array(seqs)
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    n_unique = len(set(seqs))
    print(f"unique: {n_unique}")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
