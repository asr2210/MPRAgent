"""Exp 027: random seed 27182 (e)."""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
rng = np.random.default_rng(27182)
bases = np.array(list("ACGT"))
with open(OUT, "w") as f:
    for _ in range(50_000):
        f.write("".join(bases[rng.integers(0, 4, 200)]) + "\n")
print("Done")
