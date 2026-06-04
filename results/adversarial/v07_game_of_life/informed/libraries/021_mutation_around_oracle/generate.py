"""Experiment 021 — Mutation-around-oracle library.

Take 5000 top oracle-selected seeds, generate 9 single-base-mutated variants
of each (different position for each variant) → 50k including originals.

Tests whether localized exploration of active sequence neighborhoods
provides equivalent training info to fully-diverse oracle-selected
sequences.
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
N_SEEDS = 5_000
N_VARIANTS = 10                   # 1 original + 9 mutants per seed
N_TARGET = N_SEEDS * N_VARIANTS

SEEDS_FILE = os.path.join(REPO, "libraries", "007_malinois_top_magnitude",
                          "sequences_0.txt")
BASES = "ACGT"


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    print(f"Loading 007 selections from {SEEDS_FILE}...")
    with open(SEEDS_FILE) as f:
        all_007 = [line.strip() for line in f if line.strip()]
    # Take first N_SEEDS — they were originally sorted by descending Malinois
    # score then shuffled before writing, so first N is effectively random
    # within the 007 selection. That's fine; we want a representative subset
    # of oracle-actives.
    seeds = all_007[:N_SEEDS]
    print(f"Got {len(seeds)} oracle seed sequences")

    out_seqs = []
    seen = set()
    for s in seeds:
        # Always include original
        out_seqs.append(s)
        seen.add(s)
        # Now N_VARIANTS - 1 mutants
        for _ in range(N_VARIANTS - 1):
            for _try in range(10):
                pos = int(rng.integers(0, SEQLEN))
                orig = s[pos]
                cand_bases = [b for b in BASES if b != orig]
                new_base = str(rng.choice(cand_bases))
                mutant = s[:pos] + new_base + s[pos+1:]
                if mutant not in seen:
                    out_seqs.append(mutant)
                    seen.add(mutant)
                    break

    print(f"Total: {len(out_seqs):,}")
    assert len(out_seqs) == N_TARGET, f"got {len(out_seqs)}, want {N_TARGET}"

    rng.shuffle(out_seqs)

    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in out_seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")

    # Diversity check: how many unique?
    print(f"Unique sequences: {len(set(out_seqs)):,}")

    write_sequences(out_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
