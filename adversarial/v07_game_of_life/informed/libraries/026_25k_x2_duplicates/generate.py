"""Experiment 026 — 25k unique random × 2 exact duplicates = 50k.

Test whether duplication penalty scales with redundancy depth or saturates.
- 023 (25k unique + 25k 1-mut): mean_r 0.380
- 026 (25k unique + 25k EXACT dupes): predicts somewhere 0.35-0.38

Compares same unique-count (25k) but redundancy = mutants vs exact copies.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)

from utils.seqlib import write_sequences, SEQLEN

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N_UNIQUE = 25_000
N_COPIES = 2


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    print(f"Generating {N_UNIQUE:,} unique GC-50 random sequences...")
    lookup = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_UNIQUE, SEQLEN), dtype=np.int8)
    unique = ["".join(row) for row in lookup[idx]]
    assert len(set(unique)) == N_UNIQUE

    all_seqs = unique * N_COPIES
    print(f"Total after {N_COPIES}x duplication: {len(all_seqs):,}")
    assert len(all_seqs) == 50_000
    assert len(set(all_seqs)) == N_UNIQUE

    rng.shuffle(all_seqs)
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")
    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
