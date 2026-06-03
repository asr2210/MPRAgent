#!/usr/bin/env python3
"""
Experiment 005 — Dirichlet(0.1) very-extreme composition.

Tests monotonicity of alpha in Dirichlet → eval_01 from exp 004 (alpha=0.3,
eval_01=0.0786) vs exp's-of-strategies-baseline Dirichlet(1.0) (eval_01=0.0768).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
ALPHA = 0.1
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for i in range(N_TOTAL):
        p = rng.dirichlet([ALPHA] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        seqs.append("".join(ALPHABET[idx]))
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences (Dirichlet alpha={ALPHA}) to {OUT}")


if __name__ == "__main__":
    main()
