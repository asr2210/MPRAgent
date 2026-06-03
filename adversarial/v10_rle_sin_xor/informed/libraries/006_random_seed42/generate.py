"""
Experiment 006: pure random uniform, seed 42.

Test seed-to-seed variance of mean_r. My synth_random (seed 1) got 0.5210.
Baseline random_uniform across 5 seeds was 0.5202. Different seed checks
whether single-seed numbers are reliable indicators.
"""
import os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")
rng = np.random.default_rng(42)
bases = np.array(list("ACGT"))
with open(OUT, "w") as f:
    for _ in range(50_000):
        f.write("".join(bases[rng.integers(0, 4, 200)]) + "\n")
