"""
Experiment 003: pure synthetic random sequences.

Goal: diagnostic. If synthetic random gives nonzero SK-N-SH r (as baselines
suggest synth_oracle does), then the SK-N-SH≈0 issue in experiments 001/002
is due to my DHS approach. If synthetic also gives SK-N-SH≈0, the issue is
upstream of sequence content.

Baseline: synth_oracle gets eval_01=0.6840 on 50k → expected mean r around
that level across cells.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
N_SEQ = 50_000
WINDOW = 200
SEED = 1

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))
with open(OUT, "w") as f:
    for _ in range(N_SEQ):
        f.write("".join(bases[rng.integers(0, 4, WINDOW)]) + "\n")
print("done")
