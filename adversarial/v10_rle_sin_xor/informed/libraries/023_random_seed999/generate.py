"""Exp 023: random seed 999 (seed sweep)."""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
rng = np.random.default_rng(999)
bases = np.array(list("ACGT"))
with open(OUT, "w") as f:
    for _ in range(50_000):
        f.write("".join(bases[rng.integers(0, 4, 200)]) + "\n")
print("Done")
