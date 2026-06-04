"""Experiment 023 — Independence threshold test: 25k unique random + 1 mutant each.

Tests the bridge between 020 (25k unique + 25k RC = 50k effective indep) and
021 (5k unique + 9 mutants = 5k effective indep). At 25k unique + 25k 1-base
mutants, effective independence should still be near 25k (mutants share 199/200
positions with original).

Prediction (theory v15): if effective indep ≈ 25k, mean_r ≈ 0.386 (matches 020).
If even shorter correlation kills info, will be below 0.386.
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
BASES = "ACGT"


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    print(f"Generating {N_UNIQUE:,} unique GC-50 random sequences...")
    lookup = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_UNIQUE, SEQLEN), dtype=np.int8)
    seqs = ["".join(row) for row in lookup[idx]]
    seen = set(seqs)

    print(f"Generating 1 mutant per seed...")
    mutants = []
    for s in seqs:
        for _try in range(10):
            pos = int(rng.integers(0, SEQLEN))
            orig = s[pos]
            cand = [b for b in BASES if b != orig]
            new_base = str(rng.choice(cand))
            mut = s[:pos] + new_base + s[pos+1:]
            if mut not in seen:
                mutants.append(mut)
                seen.add(mut)
                break

    all_seqs = seqs + mutants
    print(f"Total: {len(all_seqs):,}")
    assert len(all_seqs) == 50_000

    rng.shuffle(all_seqs)
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")
    print(f"Unique: {len(set(all_seqs)):,}")
    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
