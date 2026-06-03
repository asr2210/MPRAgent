#!/usr/bin/env python3
"""
Experiment 023 — Fine Dirichlet alpha mix around peak.

16,667 sequences each at alpha ∈ {0.2, 0.3, 0.4} = 50,001 → 50,000.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
ALPHAS = [0.2, 0.3, 0.4]
N_PER_ALPHA = N_TOTAL // len(ALPHAS)  # 16666
SEQ_LEN = 200
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for alpha in ALPHAS:
        n = N_PER_ALPHA + (1 if alpha == ALPHAS[0] else 0)
        if alpha == ALPHAS[0]:
            n = N_TOTAL - N_PER_ALPHA * (len(ALPHAS) - 1)
        else:
            n = N_PER_ALPHA
        # ensure total = 50000
        for _ in range(n):
            p = rng.dirichlet([alpha] * 4)
            idx = rng.choice(4, size=SEQ_LEN, p=p)
            seqs.append("".join(ALPHABET[idx]))
    # Trim/pad to exactly 50000
    if len(seqs) < N_TOTAL:
        for _ in range(N_TOTAL - len(seqs)):
            p = rng.dirichlet([0.3] * 4)
            idx = rng.choice(4, size=SEQ_LEN, p=p)
            seqs.append("".join(ALPHABET[idx]))
    seqs = seqs[:N_TOTAL]
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
