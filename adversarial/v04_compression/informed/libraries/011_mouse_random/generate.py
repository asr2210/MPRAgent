#!/usr/bin/env python3
"""
011_mouse_random — 50k random 200bp tiles from mouse mm10 chromosomes
(chr1, chr11, chr17, chr19 — gene-rich, comparable in size/density to
the human chr1/17/19/22 used in E2).

Diagnostic for cross-species. Tests whether the real-DNA prior is
vertebrate-universal (mouse should score near 0.45-0.49) or human-
specific (mouse << 0.40, similar to motif-implant or worse).

Three possible outcomes:
- mouse ~ 0.45-0.49 (human - 0.02 to -0.05):
    real-DNA prior is mostly vertebrate-universal. Cross-species
    mixing becomes a promising direction for breaking 0.50.
- mouse ~ 0.30-0.40:
    real-DNA prior is partly species-specific. Mouse is not a free
    substitute for human, but it might still help as diversity if
    mixed. E12 would test human+mouse.
- mouse < 0.30:
    real-DNA prior is heavily human-specific (likely TFBS sequences
    that diverged). Drop cross-species; pivot to oracle-stratified
    or cCRE-mix designs.

Seed 0. Length-weighted random sampling across the 4 chromosomes.
N-rejection at the window level.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
CHROMS = ["chr1", "chr11", "chr17", "chr19"]
DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "mm10"


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    print("loading mouse chromosomes...")
    chrom_seqs = {c: load_chrom(DATA_DIR / f"{c}.fa.gz") for c in CHROMS}
    lengths = {c: len(s) for c, s in chrom_seqs.items()}
    print(f"  loaded {len(chrom_seqs)} chroms, total {sum(lengths.values()):,} bp")
    for c, l in lengths.items():
        print(f"    {c}: {l:,}")

    # Length-weighted chromosome choice
    chroms_list = list(CHROMS)
    weights = np.array([lengths[c] for c in chroms_list], dtype=float)
    weights /= weights.sum()

    valid = set("ACGT")
    out = []
    attempts = 0
    while len(out) < N and attempts < N * 10:
        attempts += 1
        c = chroms_list[rng.choice(len(chroms_list), p=weights)]
        seq = chrom_seqs[c]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start : start + L]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    print(f"wrote {len(out)} of {attempts} attempts "
          f"(reject {1 - len(out)/attempts:.3f})")
    assert len(out) == N
    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
