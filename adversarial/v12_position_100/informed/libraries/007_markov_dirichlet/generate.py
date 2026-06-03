#!/usr/bin/env python3
"""
Experiment 007 — per-sequence first-order Markov chain over bases.

For each of 50,000 sequences: draw a 4x4 transition matrix (each row ~
Dirichlet(0.3) over the 4 bases), then generate 200bp via the chain from a
random starting base.

Hypothesis: this gives every sequence a unique DINUCLEOTIDE profile, on top
of the per-sequence base composition variance already captured by exp 004
(Dirichlet(0.3) at 0.0786). If dinucleotide axes are also learnable at 50k,
eval_01 should exceed 0.0786.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for _ in range(N_TOTAL):
        T = np.empty((4, 4))
        for i in range(4):
            T[i] = rng.dirichlet([ALPHA] * 4)
        idx = np.empty(SEQ_LEN, dtype=np.int8)
        # Random start
        idx[0] = rng.integers(0, 4)
        for j in range(1, SEQ_LEN):
            idx[j] = rng.choice(4, p=T[idx[j - 1]])
        seqs.append("".join(ALPHABET[idx]))

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
