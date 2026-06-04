#!/usr/bin/env python3
"""
Experiment 022 — Dirichlet(0.3) with reverse-complement augmentation.

25k unique Dirichlet(0.3) sequences + 25k of their reverse complements = 50k total.
Tests if strand-symmetry helps. Theory v5 predicts no (composition unchanged by RC
counting-wise but with permuted base labels — same training info).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_UNIQUE = 25_000
N_TOTAL = 50_000
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))
RC = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rev_complement(s):
    return "".join(RC[c] for c in reversed(s))


def main():
    rng = np.random.default_rng(SEED)
    fwd = []
    for _ in range(N_UNIQUE):
        p = rng.dirichlet([ALPHA] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        fwd.append("".join(ALPHABET[idx]))
    rev = [rev_complement(s) for s in fwd]
    all_seqs = np.array(fwd + rev)
    rng.shuffle(all_seqs)

    assert len(all_seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in all_seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(all_seqs))} (some self-RC pairs may collide)")
    with open(OUT, "w") as fh:
        for s in all_seqs:
            fh.write(s + "\n")
    print(f"wrote {len(all_seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
