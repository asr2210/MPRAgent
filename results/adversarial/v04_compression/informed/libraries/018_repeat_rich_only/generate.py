#!/usr/bin/env python3
"""
018_repeat_rich_only — INVERSE of E17. Keep only tiles that FAIL the
quality filter (repeat-rich / low-complexity).

A tile is KEPT if ANY of:
- max single-nucleotide frequency > 0.45
- 2-mer Shannon entropy < 3.0 bits
- max homopolymer run ≥ 5

Tests the hypothesis from E17: if repeat-rich is informative, this
should give a non-trivial score. If repeat-rich alone is too narrow,
this lands well below E2.

E17 (high-complexity only): 0.407 — much worse than E2 (0.499).
If E18 ≈ 0.40: complexity range needs BOTH halves, neither alone works.
If E18 > 0.45: repeat-rich IS the dominant signal.
If E18 ≈ 0.30: repeat-rich is mostly noise.
"""
import gzip
import math
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
CHROMS = ["chr1", "chr17", "chr19", "chr22"]
DATA = Path(__file__).resolve().parents[2] / "data"


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def max_nuc_freq(s: str) -> float:
    n = len(s)
    return max(s.count(b) for b in "ACGT") / n


def kmer_entropy(s: str, k: int) -> float:
    counts = {}
    for i in range(len(s) - k + 1):
        kmer = s[i:i+k]
        counts[kmer] = counts.get(kmer, 0) + 1
    total = sum(counts.values())
    h = 0.0
    for c in counts.values():
        p = c / total
        h -= p * math.log2(p)
    return h


def max_homopolymer(s: str) -> int:
    best = cur = 1
    prev = s[0]
    for c in s[1:]:
        if c == prev:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
            prev = c
    return best


def is_repeat_rich(s: str) -> bool:
    """Returns True if s is repeat-rich / low-complexity."""
    if max_nuc_freq(s) > 0.45:
        return True
    if kmer_entropy(s, 2) < 3.0:
        return True
    if max_homopolymer(s) >= 5:
        return True
    return False


def main():
    rng = np.random.default_rng(SEED)
    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}
    weights = np.array([len(chrom_seqs[c]) - L for c in CHROMS], dtype=float)
    weights /= weights.sum()

    valid = set("ACGT")
    out = []
    n_attempts = 0
    while len(out) < N and n_attempts < N * 10:
        n_attempts += 1
        ci = rng.choice(len(CHROMS), p=weights)
        seq = chrom_seqs[CHROMS[ci]]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start:start + L]
        if not set(w).issubset(valid):
            continue
        if not is_repeat_rich(w):
            continue
        out.append(w)

    print(f"attempts: {n_attempts:,}, accepted: {len(out):,}")
    assert len(out) == N

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
