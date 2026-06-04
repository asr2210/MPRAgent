"""
Experiment 005: Dirichlet with per-sequence alpha drawn from a wide range.

Pure dirichlet_composition (exp 002) used alpha = 0.5 for all sequences.
This experiment varies alpha per sequence to create a wider range of
composition concentrations:
  - small alpha (~0.3) → strongly skewed compositions
  - medium alpha (~0.5) → moderate diversity
  - large alpha (~2.0) → near-uniform compositions
Avoids alpha < 0.3 because that produces near-homopolymers
(v06 baseline homopolymer_rich = 0.058, very bad).

Hypothesis: a wider range of composition concentrations exposes the
model to more diverse k-mer-frequency profiles than a single fixed
alpha, potentially improving on dirichlet(0.5) = 0.1395.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    # log-uniform alpha between 0.3 and 2.0
    log_alpha = rng.uniform(np.log(0.3), np.log(2.0), size=N)
    alphas = np.exp(log_alpha)

    seqs = []
    for i in range(N):
        a = alphas[i]
        p = rng.dirichlet((a, a, a, a))
        b = rng.choice(4, size=LEN, p=p)
        seqs.append("".join(BASES[b]))

    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} seqs (alpha log-uniform 0.3..2.0) to {OUT}")


if __name__ == "__main__":
    main()
