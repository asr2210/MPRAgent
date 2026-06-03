#!/usr/bin/env python3
"""
008_genome_aug_motif — real genome tiles, each augmented with 1
strong TF motif planted at a random position.

50k 200bp sequences. Backbone: random tile from hg38 chr1/17/19/22
(matches 002). Augmentation: replace one ~6-15bp stretch in the
backbone with a randomly-chosen consensus motif from a curated
53-motif library (random strand).

Tests whether explicit motif signal *complements* natural-DNA
context. If 008 > 002 (0.499): augmentation works, motif diversity
is undersaturated in natural tiles. If ≈ 002: real DNA already
covers motif signal. If < 002: ectopic motifs disrupt useful
natural context.

Generalization argument: real DNA scaffolds preserve the natural
co-occurrence statistics and gene-context that the Markov/motif
experiments showed are crucial; the planted motif boosts exposure to
a controlled, diverse set of TFs across the 50k library. A model
that sees natural context *plus* explicit TF diversity should learn
both the universal grammar and specific motif vocabulary needed for
any cell type.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CHROMS = ["chr1", "chr17", "chr19", "chr22"]

MOTIFS = [
    "TATAAA", "CCAAT", "GGGGCGGGG", "TGACTCA", "TGACGTCA", "GGGAATTTCC",
    "CACGTG", "CAGCTG", "CAGGTG", "AGATAA", "TGTTTGT", "TGGACTTTG",
    "GTTAATCATTAA", "ATTGCGCAAT", "ATGCAAAT", "TTTATGAA", "TGTGGTT",
    "CAGTTG", "GGGCGTGGGCG", "AGAACAGAGTGTTCT", "TGTACA", "GCGCATGCGC",
    "CCATCTT", "CAGCAATT", "GTAACC", "CCGCGAGGAGGCAG", "TCAGCACCATG",
    "TTTCGCGC", "AAGTGA", "TTCCGGGAA", "ACCGGAAGT", "CTTTGT", "ATTGATTT",
    "TCAAGGTCA", "AGGTCATGGTCC", "TGCGTG", "CAGAC", "GGTGTGAA",
    "TGCCAA", "CCGCCATCTT", "TTGCCCAA", "AGGTCAAAGGTGACC", "GCCAAT",
    "GAAAATT", "CACCCAGCCT", "TGCAGTGCT", "CCTCAGGCT", "GACCAAT",
    "GAGGAA", "GGAA", "TAATTA", "AACCAC", "TCACGTGA",
]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        first = f.readline()
        assert first.startswith(">")
        return "".join(line.strip() for line in f).upper()


def rev_comp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in seq[::-1])


def main():
    rng = np.random.default_rng(SEED)
    chrom_seqs = {c: load_chrom(DATA_DIR / f"{c}.fa.gz") for c in CHROMS}
    weights = np.array([len(chrom_seqs[c]) - L for c in CHROMS], dtype=np.float64)
    weights /= weights.sum()
    valid = set("ACGT")

    out_lines = []
    attempts = 0
    while len(out_lines) < N:
        batch = max(N - len(out_lines), 1) * 2
        cidx = rng.choice(len(CHROMS), size=batch, p=weights)
        for ci in cidx:
            if len(out_lines) >= N:
                break
            attempts += 1
            full = chrom_seqs[CHROMS[ci]]
            start = rng.integers(0, len(full) - L)
            w = full[start : start + L]
            if not set(w).issubset(valid):
                continue
            # Plant one motif at a random position
            mi = rng.integers(0, len(MOTIFS))
            m = MOTIFS[mi]
            if rng.random() < 0.5:
                m = rev_comp(m)
            mlen = len(m)
            mstart = int(rng.integers(0, L - mlen + 1))
            w_aug = w[:mstart] + m + w[mstart + mlen :]
            assert len(w_aug) == L
            out_lines.append(w_aug)

    print(f"wrote {len(out_lines)} of {attempts} attempts "
          f"(reject {1 - len(out_lines)/attempts:.3f})")
    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out_lines))
        f.write("\n")


if __name__ == "__main__":
    main()
