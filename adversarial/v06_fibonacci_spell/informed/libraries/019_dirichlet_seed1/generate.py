"""Experiment 019: Dirichlet(0.5) seed=1 — seed scan."""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED, N, LEN = 1, 50_000, 200
BASES = np.array(list("ACGT"))

def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet((0.5,)*4, size=N)
    seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])]) for i in range(N)]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} dirichlet(0.5,seed={SEED}) to {OUT}")

if __name__ == "__main__":
    main()
