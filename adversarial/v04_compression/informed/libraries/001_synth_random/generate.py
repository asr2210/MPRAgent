#!/usr/bin/env python3
"""
001_synth_random — pipeline shakeout.

50,000 sequences of 200bp, drawn i.i.d. uniform over {A,C,G,T}.
Matches the published `synth_oracle` baseline (eval_01 ~ 0.684).
Purpose: verify prepare.py contract end-to-end, measure runtime, and
calibrate the floor in my own run before investing in DHS/motif data.
"""
import os
import numpy as np

SEED = 0
N = 50_000
L = 200
ALPHABET = np.array(list("ACGT"))

out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")

rng = np.random.default_rng(SEED)
idx = rng.integers(0, 4, size=(N, L), dtype=np.int8)
seqs = ALPHABET[idx]

with open(out_path, "w") as f:
    for row in seqs:
        f.write("".join(row.tolist()))
        f.write("\n")

print(f"Wrote {N} sequences of length {L} to {out_path}")
