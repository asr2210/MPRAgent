"""Experiment 002 — Real human genomic sequences tiled from hg38 chr19-22.

Tile chr19, 20, 21, 22 into non-overlapping 200bp windows; drop windows with
N or any non-ACGT character; drop windows that are >50% lowercase
(soft-masked, i.e. RepeatMasker hits); uniformly sample 150,000 windows
across the combined pool.

Comparison vs. experiment 001 (random baseline mean_r = 0.820):
- Does evolution-shaped sequence space teach the model more general
  regulatory grammar than uniform random?
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
SEQ_LEN = 200
SEED = 1
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
CHROMS = ("chr19", "chr20", "chr21", "chr22")
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")


def read_fasta_seq(path: str) -> str:
    """Read a single-record FASTA into one big string (keeps case)."""
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if chunks:
                    break
                continue
            chunks.append(line.rstrip())
    return "".join(chunks)


def tile_and_filter(seq: str) -> list[str]:
    valid = set("ACGTacgt")
    out: list[str] = []
    n_windows = len(seq) // SEQ_LEN
    for i in range(n_windows):
        w = seq[i * SEQ_LEN : (i + 1) * SEQ_LEN]
        if any(c not in valid for c in w):
            continue
        n_lower = sum(1 for c in w if c.islower())
        if n_lower > SEQ_LEN // 2:
            continue
        out.append(w.upper())
    return out


def main() -> None:
    pool: list[str] = []
    for c in CHROMS:
        seq = read_fasta_seq(os.path.join(DATA_DIR, f"{c}.fa"))
        kept = tile_and_filter(seq)
        print(f"{c}: {len(seq):,} bp -> {len(kept):,} clean windows", file=sys.stderr)
        pool.extend(kept)
    print(f"total pool: {len(pool):,} windows", file=sys.stderr)

    if len(pool) < N_SEQS:
        raise SystemExit(f"Only {len(pool):,} windows; need {N_SEQS:,}.")

    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(pool), size=N_SEQS, replace=False)
    sample = [pool[i] for i in idx]

    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in sample[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(sample))
        f.write("\n")

    print(f"Wrote {len(sample):,} sequences to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
