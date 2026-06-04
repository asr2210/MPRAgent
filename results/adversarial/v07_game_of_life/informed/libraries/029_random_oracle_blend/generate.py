"""Experiment 029 — 25k random + 25k 007 oracle blend.

Clean 2-way mix to test whether random and oracle blend constructively,
dilute, or just plateau. 012's 4-way mix at 0.380 was within plateau
noise. 029 isolates random + oracle in equal proportions.

Predictions:
- Random dilutes: 029 < 007 (~0.378)
- Oracle dilutes random: 029 ≈ 0.385
- Constructive blend: 029 > 0.386 (rare)
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
N_RAND = 25_000
N_ORACLE = 25_000

ORACLE_SRC = os.path.join(REPO, "libraries", "007_malinois_top_magnitude",
                          "sequences_0.txt")


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    print(f"Loading 007 oracle sequences...")
    with open(ORACLE_SRC) as f:
        all_007 = [line.strip() for line in f if line.strip()]
    # First 25k of 007 (007 was shuffled before writing → first N is random subset)
    oracle = all_007[:N_ORACLE]
    print(f"Got {len(oracle)} oracle sequences")
    seen = set(oracle)

    print(f"Generating {N_RAND:,} unique GC-50 random sequences (not in oracle)...")
    lookup = np.array(list("ACGT"))
    rand_seqs = []
    while len(rand_seqs) < N_RAND:
        batch = N_RAND - len(rand_seqs)
        idx = rng.integers(0, 4, size=(batch * 2, SEQLEN), dtype=np.int8)
        for row in lookup[idx]:
            s = "".join(row)
            if s in seen:
                continue
            seen.add(s)
            rand_seqs.append(s)
            if len(rand_seqs) >= N_RAND:
                break

    all_seqs = oracle + rand_seqs
    assert len(all_seqs) == 50_000
    assert len(set(all_seqs)) == 50_000

    rng.shuffle(all_seqs)
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")
    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
