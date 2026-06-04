"""
Experiment 014: random with each sequence forced to GC=50% exactly.

Pure random uniform has per-sequence GC variance: stdev = sqrt(0.5*0.5/200) ≈ 0.035.
So GC varies ~0.43–0.57 across sequences. This experiment removes that variance
by forcing every sequence to exactly 100 G+C and 100 A+T (in random positions).

Hypothesis: reducing per-sequence composition variance MAY help the model
generalize, or MAY hurt by removing useful variance. Tests how sensitive the
oracle is to per-sequence GC content.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
# Per-position: 50% AT, 50% GC; within each, random choice
with open(OUT, "w") as f:
    for _ in range(50_000):
        # 100 GC, 100 AT shuffled
        positions = np.zeros(200, dtype=int)
        gc_pos = rng.choice(200, size=100, replace=False)
        positions[gc_pos] = 1  # 1 = G/C, 0 = A/T
        seq = []
        for p in positions:
            if p == 0:
                seq.append("A" if rng.random() < 0.5 else "T")
            else:
                seq.append("G" if rng.random() < 0.5 else "C")
        f.write("".join(seq) + "\n")
print("Done")
