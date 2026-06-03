"""
Experiment 003: Markov-Dirichlet — per-sequence first-order Markov chain
with a Dirichlet-sampled transition matrix.

Extends dirichlet_composition (exp 002, 0.1395) by varying *dinucleotide*
frequencies as well as base frequencies. Each sequence has its own 4x4
transition matrix T (each row ~ Dirichlet(alpha)) and its own initial
distribution p0 ~ Dirichlet(alpha).

Hypothesis: if the weak model in prepare.py uses dinucleotide- or
short-k-mer features, exposing it to a wider range of those statistics
will give it more training signal than i.i.d. base composition variation.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
ALPHA_INIT = (0.5,) * 4
ALPHA_ROW = (0.5,) * 4
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)

    p0 = rng.dirichlet(ALPHA_INIT, size=N)                  # (N, 4)
    T = rng.dirichlet(ALPHA_ROW, size=N * 4).reshape(N, 4, 4)  # (N, 4, 4)
    cdf_T = np.cumsum(T, axis=2)                            # (N, 4, 4)

    # Sample initial state
    u0 = rng.random(N)
    cdf0 = np.cumsum(p0, axis=1)
    s = (u0[:, None] < cdf0).argmax(axis=1)                 # (N,)

    out = np.empty((N, LEN), dtype=np.uint8)
    out[:, 0] = s
    idx = np.arange(N)
    for t in range(1, LEN):
        u = rng.random(N)
        cdf_t = cdf_T[idx, s]                               # (N, 4)
        s = (u[:, None] < cdf_t).argmax(axis=1)
        out[:, t] = s

    seqs = ["".join(BASES[row]) for row in out]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} seqs to {OUT}")


if __name__ == "__main__":
    main()
