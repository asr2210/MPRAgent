#!/usr/bin/env python3
"""
001_ccre_random — Uniform sampling of ENCODE cCREs (V3, GRCh38).

Mirrors the dhs_random baseline (Pearson r ≈ 0.7089 on eval_01) but uses the
SCREEN cCRE registry (~1.06M elements across ~1,500 biosamples) rather than
the Meuleman DHS index. Goal: establish a clean ENCODE-cCRE baseline.

Each sequence is the 200bp window centered on a cCRE midpoint.
"""
import os
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

N_SEQ = 50_000
L = 200
SEED = 0

ACGT = set("ACGT")


def main():
    rng = random.Random(SEED)
    np.random.seed(SEED)

    # Load cCREs
    rows = []
    with open(CCRE_BED) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            # Skip non-standard chromosomes (we want main assembly only)
            if "_" in chrom or chrom in ("chrM",):
                continue
            rows.append((chrom, start, end))
    print(f"Loaded {len(rows)} cCREs (after chromosome filter)")

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

    # Pick a 200bp window from each chosen cCRE.
    # If cCRE >= 200bp, center on midpoint and take 200bp.
    # If cCRE < 200bp, pad equally around the cCRE.
    sequences = []
    attempts = 0
    max_attempts = N_SEQ * 5
    idx_pool = list(range(len(rows)))
    rng.shuffle(idx_pool)
    ptr = 0
    while len(sequences) < N_SEQ and attempts < max_attempts:
        attempts += 1
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
        seq = fa[chrom][s:e]
        seq = str(seq).upper()
        # Reject sequences containing N (heterochromatic / unassembled gaps).
        if len(seq) != L or not set(seq).issubset(ACGT):
            continue
        sequences.append(seq)

    print(f"Collected {len(sequences)} sequences after {attempts} attempts")
    assert len(sequences) == N_SEQ, f"Got {len(sequences)} sequences"

    with open(OUT, "w") as f:
        for s in sequences:
            f.write(s + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
