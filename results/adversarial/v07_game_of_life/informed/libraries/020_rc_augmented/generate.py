"""Experiment 020 — RC-augmented library (25k unique + 25k RCs = 50k).

Tests whether explicit reverse-complement doubling helps, matches, or hurts.
If the trainer does RC augmentation internally, this should be a no-op.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)

from utils.seqlib import write_sequences, SEQLEN, reverse_complement

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N_UNIQUE = 25_000


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    print(f"Generating {N_UNIQUE:,} unique GC-50 random sequences...")
    lookup = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_UNIQUE, SEQLEN), dtype=np.int8)
    unique_seqs = ["".join(row) for row in lookup[idx]]

    print("Generating RCs...")
    rcs = [reverse_complement(s) for s in unique_seqs]

    # Interleave so RC pairs are adjacent (might shuffle later)
    all_seqs = unique_seqs + rcs
    # Shuffle so trainer doesn't see RC pairs adjacent
    rng.shuffle(all_seqs)
    print(f"Total: {len(all_seqs)}")

    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")

    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
