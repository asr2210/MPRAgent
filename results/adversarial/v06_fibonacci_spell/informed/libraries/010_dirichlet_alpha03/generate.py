"""
Experiment 010: pure Dirichlet(0.3) — more extreme composition skew.

Dirichlet(0.5) gives 0.1395. Sobol-uniform simplex (alpha effectively 1.0)
gives 0.1371. Extreme compositions seem to help.

This experiment pushes further: alpha=0.3 → compositions even more
biased toward simplex edges. Each sequence is more likely to have one
or two dominant bases.

Risk: alpha too low → near-homopolymer sequences (we know they hurt:
v06 homopolymer_rich = 0.058).
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
ALPHA = (0.3, 0.3, 0.3, 0.3)
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet(ALPHA, size=N)
    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} dirichlet(alpha={ALPHA[0]}) seqs to {OUT}")


if __name__ == "__main__":
    main()
