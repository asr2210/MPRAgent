#!/usr/bin/env python3
"""
030_gene_density_repeat_stack — final stacking attempt.

Combine the two strongest confirmed signals into one design:
  (a) gene-density² weighting at 500kb bins (E28 = 0.5086)
  (b) repeat-rich content filter (E18 = 0.4937 vs E17 = 0.4069 — repeat
      content IS informative; high-complexity filter is catastrophic)

A tile is REJECTED if it is "high complexity" — i.e., max_nuc_freq ≤ 0.45
AND 2mer entropy ≥ 3.0 AND max homopolymer ≤ 4. Keeps the repeat-rich
~70% within gene-dense bins.

E28 (gene-density² only):           0.5086
E18 (repeat-rich only, chr1/17/19/22): 0.4937
This (gene-density² + repeat-rich): ?

If > 0.5086: signals stack — repeat content matters within gene-dense
  regions too, on top of bin selection.
If ≈ 0.5086: gene-density already captures the relevant repeat-rich
  bias (e.g., gene-dense regions are GC-rich and naturally pass).
If < 0.5086: filtering away high-complexity hurts when bins are already
  gene-biased.

3-seed.
"""
import gzip
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

N = 50_000
L = 200
BIN_SIZE = 500_000
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
DATA = Path(__file__).resolve().parents[2] / "data"
EPS = 0.5
EXP = 2.0
SEEDS = [0, 1, 2]


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
    if max_nuc_freq(s) > 0.45:
        return True
    if kmer_entropy(s, 2) < 3.0:
        return True
    if max_homopolymer(s) >= 5:
        return True
    return False


def main():
    seen = set()
    gene_starts = defaultdict(list)
    with gzip.open(DATA / "hg38.refGene.gtf.gz", "rt") as f:
        for line in f:
            parts = line.split("\t")
            if parts[2] != "transcript":
                continue
            chrom = parts[0]
            if chrom not in set(CHROMS):
                continue
            t = int(parts[3]) if parts[6] == "+" else int(parts[4])
            key = (chrom, t)
            if key in seen:
                continue
            seen.add(key)
            gene_starts[chrom].append(t)

    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in CHROMS}

    bin_weights = []
    for c in CHROMS:
        chrom_len = len(chrom_seqs[c])
        n_bins = chrom_len // BIN_SIZE
        counts = np.zeros(n_bins, dtype=np.int32)
        for g in gene_starts.get(c, []):
            b = g // BIN_SIZE
            if 0 <= b < n_bins:
                counts[b] += 1
        for b in range(n_bins):
            bin_weights.append((c, b, (counts[b] + EPS) ** EXP))

    weights = np.array([w for _, _, w in bin_weights], dtype=float)
    weights /= weights.sum()
    print(f"{len(bin_weights):,} bins; top-50 mass: {sum(np.sort(weights)[-50:]):.4f}")

    valid = set("ACGT")
    out_dir = Path(__file__).resolve().parent
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        out = []
        attempts = 0
        while len(out) < N and attempts < N * 20:
            attempts += 1
            bi = rng.choice(len(bin_weights), p=weights)
            c, b, _ = bin_weights[bi]
            seq = chrom_seqs[c]
            bin_start = b * BIN_SIZE
            bin_end = min(bin_start + BIN_SIZE, len(seq))
            if bin_end - bin_start < L:
                continue
            start = int(rng.integers(bin_start, bin_end - L))
            w = seq[start:start + L]
            if not set(w).issubset(valid):
                continue
            if not is_repeat_rich(w):
                continue
            out.append(w)
        assert len(out) == N, f"only got {len(out)} after {attempts} attempts"
        out_path = out_dir / f"sequences_{seed}.txt"
        with open(out_path, "w") as f:
            f.write("\n".join(out))
            f.write("\n")
        print(f"seed {seed}: wrote {len(out)} in {attempts:,} attempts")


if __name__ == "__main__":
    main()
