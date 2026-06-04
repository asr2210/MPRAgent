"""
Experiment 007: Sobol quasi-random simplex coverage.

Replaces independent Dirichlet draws with a Sobol low-discrepancy
sequence mapped to the 4-simplex (via exponential trick), giving more
uniform coverage of the composition space across 50,000 sequences.

If composition simplex coverage is the bottleneck, this should beat
dirichlet(0.5) = 0.1395. If dirichlet already covers well, this will
be a wash.
"""
from pathlib import Path

import numpy as np
from scipy.stats import qmc

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)

    sobol = qmc.Sobol(d=4, scramble=True, seed=SEED)
    u = sobol.random(N)
    u = np.clip(u, 1e-10, 1.0)
    e = -np.log(u)
    probs = e / e.sum(axis=1, keepdims=True)

    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))

    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} seqs (Sobol simplex) to {OUT}")


if __name__ == "__main__":
    main()
