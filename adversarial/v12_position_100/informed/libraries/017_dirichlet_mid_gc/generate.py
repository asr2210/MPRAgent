#!/usr/bin/env python3
"""
Experiment 017 — Dirichlet(0.3) with mid-GC rejection.

Sample Dirichlet(0.3) compositions, keep only those with GC ∈ [0.3, 0.7].
Tests if mid-GC concentration matches eval target distribution better than
the natural Dirichlet(0.3) which extends to extreme GC.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
ALPHA = 0.3
GC_LOW, GC_HIGH = 0.3, 0.7
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    accepted = 0
    tried = 0
    while accepted < N_TOTAL:
        # batch sample for speed
        batch = rng.dirichlet([ALPHA] * 4, size=2000)  # order A,C,G,T
        gc = batch[:, 1] + batch[:, 2]
        mask = (gc >= GC_LOW) & (gc <= GC_HIGH)
        valid = batch[mask]
        tried += 2000
        for p in valid:
            if accepted >= N_TOTAL:
                break
            idx = rng.choice(4, size=SEQ_LEN, p=p)
            seqs.append("".join(ALPHABET[idx]))
            accepted += 1
    print(f"accepted {accepted} from {tried} (acceptance rate {accepted/tried:.3f})")
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
