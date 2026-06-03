"""Experiment 016 — Bimodal per-sequence GC, library mean GC=0.50.

Half the library (25k) has per-sequence GC ≈ 0.30, the other half ≈ 0.70.
Library mean GC = 0.50, but per-sequence GC variance is large.

Tests: does the model care about per-sequence composition, or only
library-aggregate composition? If 016 ≈ 007 random GC-50, library mean
is sufficient. If 016 << 007, per-sequence composition is critical.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)

from utils.seqlib import write_sequences, SEQLEN, is_clean

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N_HALF = 25_000
GC_LOW = 0.30
GC_HIGH = 0.70


def gen_at_gc(rng, n, gc_frac):
    """Generate n sequences with per-base prob (gc_frac/2, (1-gc_frac)/2, ...)."""
    # p(C) = p(G) = gc/2; p(A) = p(T) = (1-gc)/2
    # Base order ACGT → p = [(1-gc)/2, gc/2, gc/2, (1-gc)/2]
    p = np.array([(1-gc_frac)/2, gc_frac/2, gc_frac/2, (1-gc_frac)/2])
    bases = np.array(list("ACGT"))
    # Vectorized per-sequence sampling
    idx = rng.choice(4, size=(n, SEQLEN), p=p)
    return ["".join(row) for row in bases[idx]]


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    print(f"Generating {N_HALF:,} sequences with GC={GC_LOW}...")
    low_seqs = gen_at_gc(rng, N_HALF, GC_LOW)
    print(f"Generating {N_HALF:,} sequences with GC={GC_HIGH}...")
    high_seqs = gen_at_gc(rng, N_HALF, GC_HIGH)

    all_seqs = low_seqs + high_seqs
    rng.shuffle(all_seqs)

    # Stats
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"Library GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}, "
          f"min={gcs.min():.3f}, max={gcs.max():.3f}")
    print(f"  Half low: mean={gcs[gcs<0.5].mean():.3f}")
    print(f"  Half high: mean={gcs[gcs>=0.5].mean():.3f}")

    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
