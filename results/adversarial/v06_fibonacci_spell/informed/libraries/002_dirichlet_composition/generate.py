"""
Experiment 002: Pure dirichlet_composition (replicate v06 best baseline).

Each sequence's base frequencies are drawn from a Dirichlet(alpha) prior,
then 200 bases are sampled iid from that multinomial. Different sequences
have very different compositions (some A-rich, some GC-rich, some balanced).

v06 baseline reports eval_01 = 0.1383 for this strategy across 5 seeds.
This experiment is the single-seed replica to calibrate my own scores
against the v06 table and confirm reproducibility.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
ALPHA = (0.5, 0.5, 0.5, 0.5)  # symmetric Dirichlet, spread compositions
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet(ALPHA, size=N)
    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} seqs to {OUT}")


if __name__ == "__main__":
    main()
