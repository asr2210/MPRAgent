"""
Experiment 012: pure random uniform, seed=7. Extra variance datapoint.

Current best: 0.5221 (seed 42). Seed 1: 0.5210. Seed variance ~0.001.
Run a 3rd seed to bound the variance distribution and identify the "ceiling"
random uniform actually achieves with any seed.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 7

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))
with open(OUT, "w") as f:
    for _ in range(50_000):
        f.write("".join(bases[rng.integers(0, 4, 200)]) + "\n")
print("Done")
