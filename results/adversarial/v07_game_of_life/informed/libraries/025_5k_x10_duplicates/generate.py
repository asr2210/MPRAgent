"""Experiment 025 — 5k unique random × 10 exact duplicates = 50k.

Isolates (A) effective-independence vs (B) spatial-correlation as drivers
of 021's 0.026 drop. 021 uses mutants (spatially correlated). 025 uses
exact duplicates (zero-correlation, just redundant).

Predictions (v15):
- If indep alone drives gap: 025 ≈ 0.360 (like 021)
- If spatial correlation is the mechanism: 025 ≈ 0.380+ (closer to plateau)
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
N_UNIQUE = 5_000
N_COPIES = 10


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    print(f"Generating {N_UNIQUE:,} unique GC-50 random sequences...")
    lookup = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_UNIQUE, SEQLEN), dtype=np.int8)
    unique = ["".join(row) for row in lookup[idx]]
    print(f"Unique: {len(set(unique)):,}")
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
