#!/usr/bin/env python3
"""
Experiment 027 — Multi-seed Dirichlet(0.3) ensemble.

~16,667 sequences each at seeds 0, 1, 2 = ~50k. Tests if multi-seed sampling
gives better effective composition coverage than single-seed.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
N_TOTAL = 50_000
SEEDS = [0, 1, 2]
N_PER_SEED = N_TOTAL // len(SEEDS)  # 16666
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    seqs = []
    for i, s in enumerate(SEEDS):
        rng = np.random.default_rng(s)
        n = N_PER_SEED + (N_TOTAL - N_PER_SEED * len(SEEDS) if i == 0 else 0)
        for _ in range(n):
            p = rng.dirichlet([ALPHA] * 4)
            idx = rng.choice(4, size=SEQ_LEN, p=p)
            seqs.append("".join(ALPHABET[idx]))
    final_rng = np.random.default_rng(99)
    arr = np.array(seqs)
    final_rng.shuffle(arr)

    assert len(arr) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in arr if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(arr))}")
    with open(OUT, "w") as fh:
        for s in arr:
            fh.write(s + "\n")
    print(f"wrote {len(arr)} sequences to {OUT}")


if __name__ == "__main__":
    main()
