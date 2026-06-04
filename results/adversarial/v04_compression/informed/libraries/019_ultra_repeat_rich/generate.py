#!/usr/bin/env python3
"""
019_ultra_repeat_rich — even more extreme repeat content than E18.

Keep only tiles satisfying AT LEAST ONE of:
- max single-nucleotide frequency > 0.55 (more biased than E18's 0.45)
- 2-mer Shannon entropy < 2.5 (much lower than E18's 3.0)
- max homopolymer run ≥ 8 (stricter than E18's 5)

Tests whether the "repeat-rich helps" finding from E18 is monotonic.

E2 random:               0.4992
E18 repeat-rich (70%):   0.4937
E17 high-complexity:     0.4069

Outcomes:
- E19 > E18: monotonic; even more repeat-rich = even better. Push further.
- E19 ≈ E18: saturated within repeat-rich.
- E19 < E18: diminishing returns; the moderate-repeat tier in E18 was
  doing some work.
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


def is_ultra_repeat(s: str) -> bool:
    if max_nuc_freq(s) > 0.55:
        return True
    if kmer_entropy(s, 2) < 2.5:
        return True
    if max_homopolymer(s) >= 8:
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
    while len(out) < N and n_attempts < N * 100:
        n_attempts += 1
        ci = rng.choice(len(CHROMS), p=weights)
        seq = chrom_seqs[CHROMS[ci]]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start:start + L]
        if not set(w).issubset(valid):
            continue
        if not is_ultra_repeat(w):
            continue
        out.append(w)

    print(f"attempts: {n_attempts:,}, accepted: {len(out):,} "
          f"(rate {len(out)/n_attempts:.3f})")
    assert len(out) == N

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
