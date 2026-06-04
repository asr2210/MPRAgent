"""Experiment 027 — RC-augmented oracle: 25k oracle-top + their 25k RCs.

Mirror of 020 (RC-augmented random) but with oracle-selected seeds.
Tests whether RC-equivariance holds for ACTIVE (motif-rich, often
strand-asymmetric) sequences, not just random.

If trainer's RC-equivariance is universal: matches 007 (~0.386).
If active motifs break RC-equivariance: 027 > 007 (free info via RC).
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
N_BASE = 25_000

SOURCE = os.path.join(REPO, "libraries", "007_malinois_top_magnitude",
                     "sequences_0.txt")


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    print(f"Loading 007 oracle seeds from {SOURCE}...")
    with open(SOURCE) as f:
        all_007 = [line.strip() for line in f if line.strip()]
    print(f"Got {len(all_007)} oracle sequences")
    # Take first 25k (007 was shuffled before writing → first N is effectively random subset)
    seeds = all_007[:N_BASE]
    print(f"Selected first {len(seeds)} for RC augmentation")

    print("Generating RCs...")
    rcs = [reverse_complement(s) for s in seeds]

    # Check for collisions (extremely rare for 200-mers)
    union_set = set(seeds) | set(rcs)
    print(f"Unique union: {len(union_set):,} (out of {2*N_BASE} =, expected 50000)")

    all_seqs = seeds + rcs
    assert len(all_seqs) == 50_000

    rng.shuffle(all_seqs)
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")
    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
