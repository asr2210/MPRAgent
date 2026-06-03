"""Experiment 023: Dirichlet(0.55) — fine alpha sweep around 0.5."""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED, N, LEN, ALPHA = 42, 50_000, 200, 0.55
BASES = np.array(list("ACGT"))

def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet((ALPHA,)*4, size=N)
    seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])]) for i in range(N)]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} dirichlet(alpha={ALPHA}) to {OUT}")

if __name__ == "__main__":
    main()
