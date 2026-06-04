"""
Experiment 001: Random uniform DNA calibration baseline.

50,000 sequences of 200bp, each base i.i.d. uniform {A,C,G,T}.
Equivalent to the baseline `synth_oracle` strategy. Single seed for
single-shot pipeline validation and seed-variance estimation.

Expected eval_01 ≈ 0.6840 (5-seed mean from instructions.md).
"""

import numpy as np
import os

SEED = 1
N_SEQS = 50_000
LEN = 200
BASES = np.array(["A", "C", "G", "T"])

rng = np.random.default_rng(SEED)
idx = rng.integers(0, 4, size=(N_SEQS, LEN))
seqs = np.array(["".join(BASES[row]) for row in idx])

out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out_path, "w") as f:
    f.write("\n".join(seqs.tolist()))
    f.write("\n")

print(f"Wrote {N_SEQS} sequences of length {LEN} to {out_path}")
