"""
Experiment 018: extreme bimodal GC: half GC=0.2, half GC=0.8.

More extreme version of 017. Tests if more compositional separation helps.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
with open(OUT, "w") as f:
    for i in range(50_000):
        gc_target = 0.2 if i < 25_000 else 0.8
        is_gc = rng.random(200) < gc_target
        seq = []
        for ig in is_gc:
            if ig:
                seq.append("G" if rng.random() < 0.5 else "C")
            else:
                seq.append("A" if rng.random() < 0.5 else "T")
        f.write("".join(seq) + "\n")
print("Done")
