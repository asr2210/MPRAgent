"""
Experiment 024: Mixture of Dirichlet(0.5) backbone + Dirichlet(0.2) extreme tail.

45k from Dirichlet(0.5) + 5k from Dirichlet(0.2). Tests whether a small
dose of EXTREME compositions (near-homopolymer) added to the optimal
Dirichlet(0.5) library boosts K562 (which prefers extreme compositions,
exp 013 confirmed).

Pure Dirichlet(0.2) compositions are very skewed (max-base ~0.75). At
5%, they shouldn't dominate the library but might lift K562 head.
"""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
N_MAIN = 45_000
N_EXTREME = 5_000
LEN = 200
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs_main = rng.dirichlet((0.5,)*4, size=N_MAIN)
    probs_extreme = rng.dirichlet((0.2,)*4, size=N_EXTREME)
    probs = np.concatenate([probs_main, probs_extreme])
    rng.shuffle(probs)
    seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])]) for i in range(N)]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} ({N_MAIN} alpha=0.5 + {N_EXTREME} alpha=0.2) to {OUT}")


if __name__ == "__main__":
    main()
