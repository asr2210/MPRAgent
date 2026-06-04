#!/usr/bin/env python3
"""
Experiment 009 — Replicated Dirichlet for noise averaging test.

5,000 unique Dirichlet(0.3) sequences × 10 copies each = 50,000 lines.

If the prepare.py harness measures each line independently and the trained
model averages duplicate (x,y) pairs, this reduces effective per-sequence
noise by ~sqrt(10) ≈ 3×. Costs: 90% less unique-sequence coverage.

Outcomes:
- Big win (>0.10 on eval_01): noise is the bottleneck; replication is a lever.
- Same (~0.078): noise reduction doesn't compensate for coverage loss.
- Big loss (<0.05): harness dedupes or sequence coverage dominates.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_UNIQUE = 5_000
N_COPIES = 10
N_TOTAL = N_UNIQUE * N_COPIES
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))

assert N_TOTAL == 50_000


def main():
    rng = np.random.default_rng(SEED)
    uniques = []
    for _ in range(N_UNIQUE):
        p = rng.dirichlet([ALPHA] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        uniques.append("".join(ALPHABET[idx]))

    # Replicate and shuffle
    seqs = np.array(uniques * N_COPIES, dtype=object)
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"

    n_unique_actual = len(set(seqs))
    print(f"unique sequences in file: {n_unique_actual} (intended {N_UNIQUE})")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
