"""
Experiment 014: Dirichlet(0.5) with seed=99 (vs seed=42 for exp 002).

Pure replica of exp 002 (best baseline = 0.1395) but with a different
random seed. Purpose: measure single-seed noise floor.

If 014 in [0.135, 0.145] → noise is ~0.005, 002's 0.1395 is stable
If 014 > 0.140 or < 0.135 → seed=42 was lucky/unlucky; noise is larger

Critical for interpreting all my deviations from 0.1395 baseline. Most
experiments fall in [0.1314, 0.1395]; knowing the noise floor lets me
distinguish "real losses" from "sampling noise".
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 99
N = 50_000
LEN = 200
ALPHA = (0.5, 0.5, 0.5, 0.5)
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet(ALPHA, size=N)
    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} dirichlet(alpha=0.5,seed=99) seqs to {OUT}")


if __name__ == "__main__":
    main()
