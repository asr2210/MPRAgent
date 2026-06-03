"""
Experiment 015: random with each sequence forced to exact 50/50/50/50 ACGT.

Even stricter than 014 — every sequence has exactly 50 of each base.
Removes ALL per-sequence compositional variance. Maximally "balanced" random.

Hypothesis: stripping per-sequence composition variance should either be
neutral (if oracle doesn't care about composition) or hurt (less diverse input).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
template = np.array(list("A" * 50 + "C" * 50 + "G" * 50 + "T" * 50))
with open(OUT, "w") as f:
    for _ in range(50_000):
        perm = rng.permutation(200)
        f.write("".join(template[perm]) + "\n")
print("Done")
