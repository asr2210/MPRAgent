"""Experiment 013 — Pure RC-augmented cCRE (no synthetic component).

Isolates the RC contribution from the synthetic contribution by removing
the 15k motif-embedded synthetic block from exp 012 and replacing those
slots with another 7.5k unique cCREs (also × 2 strands).

Composition:
- 75,000 unique cCRE-centered 200bp windows (forward), seed=12
- 75,000 reverse complements of those 75,000 cCREs

If pure RC reaches ≈ 0.888 (where pure cCRE + synthetic was), the synthetic
contribution is redundant once we have RC. If pure RC is significantly
below 012 (0.890), RC and synthetic contribute additively.
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_UNIQUE = 75_000
N_SEQS = 2 * N_UNIQUE  # 150,000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 12
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")

RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def read_fasta_seq(path: str) -> str:
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if chunks:
                    break
                continue
            chunks.append(line.rstrip())
    return "".join(chunks)


def build_ccre_pool() -> list[str]:
    chroms_needed: set[str] = set()
    cres: list[tuple[str, int]] = []
    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            mid = (int(parts[1]) + int(parts[2])) // 2
            cres.append((chrom, mid))
            chroms_needed.add(chrom)
    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        path = os.path.join(DATA_DIR, f"{c}.fa")
        if os.path.exists(path):
            chrom_seq[c] = read_fasta_seq(path)
    valid = set("ACGTacgt")
    kept: list[str] = []
    for chrom, mid in cres:
        seq = chrom_seq.get(chrom)
        if seq is None:
            continue
        lo, hi = mid - HALF, mid + HALF
        if lo < 0 or hi > len(seq):
            continue
        w = seq[lo:hi]
        if any(c not in valid for c in w):
            continue
        if sum(1 for c in w if c.islower()) > SEQ_LEN // 2:
            continue
        kept.append(w.upper())
    return kept


def main() -> None:
    rng = np.random.default_rng(SEED)
    ccre_pool = build_ccre_pool()
    print(f"cCRE pool: {len(ccre_pool):,}", file=sys.stderr)
    assert len(ccre_pool) >= N_UNIQUE
    idx = rng.choice(len(ccre_pool), size=N_UNIQUE, replace=False)
    fwd = [ccre_pool[i] for i in idx]
    rc = [revcomp(s) for s in fwd]
    combined = fwd + rc
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])
    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(f"Wrote {N_SEQS:,} ({N_UNIQUE} fwd + {N_UNIQUE} RC) to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
