"""Experiment 001 — Random baseline.

Generate 150,000 uniformly random 200bp ACGT sequences. Seed fixed so the
library is reproducible. This establishes the floor: how well can a model
learn regulatory grammar from sequences that have essentially no enrichment
for functional elements?
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
SEQ_LEN = 200
ALPHABET = np.array(list("ACGT"))
SEED = 0
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")


def main() -> None:
    rng = np.random.default_rng(SEED)
    # Draw all bases at once: shape (N_SEQS, SEQ_LEN)
    idx = rng.integers(0, 4, size=(N_SEQS, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]  # shape (N_SEQS, SEQ_LEN), dtype <U1
    # Join rows into strings without per-base Python overhead.
    seqs = chars.view(f"<U{SEQ_LEN}").ravel()
    assert seqs.shape == (N_SEQS,)
    assert all(len(s) == SEQ_LEN for s in seqs[:5])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(seqs.tolist()))
        f.write("\n")

    print(f"Wrote {N_SEQS} sequences of length {SEQ_LEN} to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
