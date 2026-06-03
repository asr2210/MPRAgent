#!/usr/bin/env python3
"""
002_ccre_replicate_x10 — Same cCREs as 001, but only 5,000 unique sequences,
each emitted 10 times (50,000 total).

Test: does replication help under v14's apparently very brief training?
If each unique pattern is seen 10x in one epoch the model gets denser per-
pattern feedback, and label noise averages within-group.
"""
import random
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
CCRE_BED = DATA / "GRCh38-cCREs.bed"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

N_UNIQUE = 5_000
REPLICATES = 10
L = 200
SEED = 0
ACGT = set("ACGT")


def main():
    rng = random.Random(SEED)
    np.random.seed(SEED)

    rows = []
    with open(CCRE_BED) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if "_" in chrom or chrom == "chrM":
                continue
            rows.append((chrom, start, end))

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

    sequences = []
    idx_pool = list(range(len(rows)))
    rng.shuffle(idx_pool)
    ptr = 0
    while len(sequences) < N_UNIQUE:
        if ptr >= len(idx_pool):
            rng.shuffle(idx_pool)
            ptr = 0
        chrom, start, end = rows[idx_pool[ptr]]
        ptr += 1
        mid = (start + end) // 2
        s = mid - L // 2
        e = s + L
        if s < 0:
            continue
        try:
            chrom_len = len(fa[chrom])
        except KeyError:
            continue
        if e > chrom_len:
            continue
        seq = str(fa[chrom][s:e]).upper()
        if len(seq) != L or not set(seq).issubset(ACGT):
            continue
        sequences.append(seq)

    # Emit each unique sequence REPLICATES times, in interleaved order so
    # replicates are spread across the training stream.
    all_seq = []
    for r in range(REPLICATES):
        for s in sequences:
            all_seq.append(s)
    assert len(all_seq) == N_UNIQUE * REPLICATES

    with open(OUT, "w") as f:
        for s in all_seq:
            f.write(s + "\n")
    print(f"Wrote {OUT}: {N_UNIQUE} unique × {REPLICATES} reps = {len(all_seq)} lines")


if __name__ == "__main__":
    main()
