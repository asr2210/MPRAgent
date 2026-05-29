"""Experiment 010 — cCRE neighbourhoods.

For each of 50k uniformly-sampled cCREs, generate 3 windows:
  center    [mid-100, mid+100)
  upstream  [mid-300, mid-100)
  downstream [mid+100, mid+300)

Apply per-window filter (ACGT only, <50% softmasked). If a window fails,
drop just that window (keep the others). To guarantee 150k sequences in
total, draw cCREs from a larger pool until exactly 150k pass.

Tests whether sequences *adjacent to* cCREs (flanks containing chromatin
context, TFBS clusters, regulatory dinucleotide composition) carry
additional generalisable signal beyond cCRE centres alone.
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 9
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")


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


def passes(w: str) -> bool:
    valid = set("ACGTacgt")
    if any(c not in valid for c in w):
        return False
    n_lower = sum(1 for c in w if c.islower())
    return n_lower <= SEQ_LEN // 2


def main() -> None:
    # Load cCREs.
    chroms_needed: set[str] = set()
    cres: list[tuple[str, int]] = []
    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            mid = (int(parts[1]) + int(parts[2])) // 2
            cres.append((chrom, mid))
            chroms_needed.add(chrom)
    print(f"Total cCREs: {len(cres):,}", file=sys.stderr)

    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        path = os.path.join(DATA_DIR, f"{c}.fa")
        if os.path.exists(path):
            chrom_seq[c] = read_fasta_seq(path)

    # Shuffle cCRE order with seed, then walk through until we have 150k windows.
    rng = np.random.default_rng(SEED)
    order = rng.permutation(len(cres))

    out: list[str] = []
    n_dropped_oob = 0
    n_dropped_filter = 0
    n_used_ccres = 0
    for ii in order:
        chrom, mid = cres[int(ii)]
        seq = chrom_seq.get(chrom)
        if seq is None:
            continue
        n_used_ccres += 1
        for offset in (0, -200, +200):
            lo = mid + offset - HALF
            hi = lo + SEQ_LEN
            if lo < 0 or hi > len(seq):
                n_dropped_oob += 1
                continue
            w = seq[lo:hi]
            if not passes(w):
                n_dropped_filter += 1
                continue
            out.append(w.upper())
            if len(out) >= N_SEQS:
                break
        if len(out) >= N_SEQS:
            break

    print(
        f"used {n_used_ccres:,} cCREs to produce {len(out):,} windows | "
        f"dropped oob={n_dropped_oob:,} filter={n_dropped_filter:,}",
        file=sys.stderr,
    )
    assert len(out) == N_SEQS, f"expected {N_SEQS}, got {len(out)}"

    rng.shuffle(out)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in out[:10])
    with open(OUT_PATH, "w") as f:
        f.write("\n".join(out))
        f.write("\n")
    print(f"Wrote {N_SEQS:,} sequences to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
