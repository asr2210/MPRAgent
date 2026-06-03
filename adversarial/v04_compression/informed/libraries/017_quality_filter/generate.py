#!/usr/bin/env python3
"""
017_quality_filter — 50k random 200bp tiles from chr1/17/19/22, but
with QUALITY FILTERING to exclude low-complexity / pure-repeat tiles.

Filter rules (reject if ANY of):
- max single-nucleotide frequency > 0.45 (catches AAAA... / TTTT...
  / simple monomer repeats)
- 2-mer Shannon entropy < 3.0 bits (out of max log2(16) = 4.0;
  catches tandem-2mer like (AC)50 = entropy 1.0)
- max 5+ consecutive identical bases (homopolymers)

Hypothesis: median random tile may include 10-30% "wasted slots" in
simple-repeat regions that provide no informative signal. Replacing
these with higher-complexity random should help.

Outcomes vs E2 (0.4992):
- ≥0.51: quality filter lifts ceiling.
- ~0.50: neutral; existing random is fine.
- <0.49: filter is too aggressive, hurts diversity.
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


def passes_filter(s: str) -> bool:
    if max_nuc_freq(s) > 0.45:
        return False
    if kmer_entropy(s, 2) < 3.0:
        return False
    if max_homopolymer(s) >= 5:
        return False
    return True


def main():
    rng = np.random.default_rng(SEED)
    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}
    weights = np.array([len(chrom_seqs[c]) - L for c in CHROMS], dtype=float)
    weights /= weights.sum()
    print(f"chrom weights: {dict(zip(CHROMS, weights))}")

    valid = set("ACGT")
    out = []
    n_attempts = 0
    n_n_reject = 0
    n_filter_reject = 0
    while len(out) < N and n_attempts < N * 50:
        n_attempts += 1
        ci = rng.choice(len(CHROMS), p=weights)
        seq = chrom_seqs[CHROMS[ci]]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start:start + L]
        if not set(w).issubset(valid):
            n_n_reject += 1
            continue
        if not passes_filter(w):
            n_filter_reject += 1
            continue
        out.append(w)

    print(f"attempts: {n_attempts:,}")
    print(f"  N-rejects: {n_n_reject:,} ({n_n_reject/n_attempts:.3f})")
    print(f"  filter-rejects: {n_filter_reject:,} ({n_filter_reject/n_attempts:.3f})")
    print(f"  accepted: {len(out):,}")
    assert len(out) == N

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
