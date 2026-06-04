"""
Experiment 028: Asymmetric Dirichlet (0.4, 0.6, 0.6, 0.4) at seed=99.

Replicates exp 026 (which scored 0.1385 at seed=42) with seed=99.
Tests if the asymmetric prior is genuinely better than symmetric or
if 026 was just seed=42 luck.

If 028 > 0.137 (mean Dirichlet noise floor): asymmetric helps regardless
of seed → real signal worth pursuing.
If 028 ~ 0.136: 026 was seed=42 luck, asymmetric is no different.
"""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 99
N = 50_000
LEN = 200
ALPHA = (0.4, 0.6, 0.6, 0.4)
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet(ALPHA, size=N)
    seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])]) for i in range(N)]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} asymmetric dirichlet seed=99 to {OUT}")


if __name__ == "__main__":
    main()
