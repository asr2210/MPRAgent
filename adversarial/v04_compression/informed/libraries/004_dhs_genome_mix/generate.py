#!/usr/bin/env python3
"""
004_dhs_genome_mix — direct test of DHS↔genome complementarity.

25k random hg38 genome tiles (chr1/17/19/22, same as 002) + 25k DHS-
stratified sequences from SynthSeqs train (same source as 003).

Decision logic for follow-ups:
- If 004 ≈ 002 (0.499): DHS adds no signal on top of broad genome.
  Pursue real-genome-broadening designs.
- If 004 > 002: DHS complements; find optimal mix proportion.
- If 004 < 002: DHS actively dilutes. Drop DHS entirely.

Generalization argument: a model trained on the mix sees both
universal genomic priors (real DNA distribution, k-mer/CpG/repeat
content) AND every NMF accessibility program equally. If DHS programs
encode regulatory grammar shared across cell types, the mix should
be at least as informative as either component alone.
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N = 50_000
N_GENOME = N // 2
N_DHS = N - N_GENOME
PER_COMPONENT = N_DHS // 16  # 1562 → 16 * 1562 = 24992; we'll fill remainder
L = 200
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CHROMS = ["chr1", "chr17", "chr19", "chr22"]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        first = f.readline()
        assert first.startswith(">")
        return "".join(line.strip() for line in f).upper()


def sample_genome(rng, n):
    chrom_seqs = {c: load_chrom(DATA_DIR / f"{c}.fa.gz") for c in CHROMS}
    weights = np.array([len(chrom_seqs[c]) - L for c in CHROMS], dtype=np.float64)
    weights /= weights.sum()
    valid = set("ACGT")
    out = []
    while len(out) < n:
        batch = max(n - len(out), 1) * 2
        cidx = rng.choice(len(CHROMS), size=batch, p=weights)
        for ci in cidx:
            if len(out) >= n:
                break
            seq_full = chrom_seqs[CHROMS[ci]]
            start = rng.integers(0, len(seq_full) - L)
            w = seq_full[start : start + L]
            if not set(w).issubset(valid):
                continue
            out.append(w)
    return out


def sample_dhs(rng, n):
    df = pd.read_csv(DATA_DIR / "synthseqs_train.csv.gz", sep="\t")
    components = sorted(df["component"].unique())
    per = n // len(components)
    rem = n - per * len(components)
    picks = []
    for i, c in enumerate(components):
        sub = df[df["component"] == c]
        take = per + (1 if i < rem else 0)
        idx = rng.choice(len(sub), size=take, replace=False)
        picks.extend(sub.iloc[idx]["raw_sequence"].tolist())
    return picks


def main():
    rng = np.random.default_rng(SEED)
    g = sample_genome(rng, N_GENOME)
    d = sample_dhs(rng, N_DHS)
    seqs = g + d
    rng.shuffle(seqs)
    assert len(seqs) == N, len(seqs)
    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s)
            f.write("\n")
    print(f"wrote {N} sequences ({N_GENOME} genome + {N_DHS} DHS)")


if __name__ == "__main__":
    main()
