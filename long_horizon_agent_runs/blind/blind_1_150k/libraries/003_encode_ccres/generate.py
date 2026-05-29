"""Experiment 003 — ENCODE V3 cCRE-centered 200bp windows.

Builds a library of 150,000 200bp sequences, each centered on the midpoint
of an ENCODE V3 candidate cis-regulatory element (cCRE). This sharply
enriches the library for putative regulatory function compared to plain
genomic tiling (experiment 002).

Comparison points:
- Experiment 001 (uniform random):   mean_r 0.820
- Experiment 002 (genomic tiles):    mean_r 0.873
- Experiment 003 (cCRE-centered):    ?
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 2
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


def main() -> None:
    # Load all hg38 chromosomes we need (lazy: only what cCREs reference).
    chroms_needed: set[str] = set()
    cres: list[tuple[str, int]] = []  # (chrom, midpoint)
    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            mid = (start + end) // 2
            cres.append((chrom, mid))
            chroms_needed.add(chrom)
    print(f"Loaded {len(cres):,} cCREs across {len(chroms_needed)} chromosomes", file=sys.stderr)

    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        path = os.path.join(DATA_DIR, f"{c}.fa")
        if not os.path.exists(path):
            print(f"WARN: missing {path} — skipping cCREs on {c}", file=sys.stderr)
            continue
        chrom_seq[c] = read_fasta_seq(path)
        print(f"  {c}: {len(chrom_seq[c]):,} bp", file=sys.stderr)

    # Build candidate windows, applying the same filter as exp 002.
    valid = set("ACGTacgt")
    kept: list[str] = []
    n_oob = n_bad_char = n_repeat = 0
    for chrom, mid in cres:
        seq = chrom_seq.get(chrom)
        if seq is None:
            continue
        lo, hi = mid - HALF, mid + HALF
        if lo < 0 or hi > len(seq):
            n_oob += 1
            continue
        w = seq[lo:hi]
        if len(w) != SEQ_LEN:
            n_oob += 1
            continue
        if any(c not in valid for c in w):
            n_bad_char += 1
            continue
        n_lower = sum(1 for c in w if c.islower())
        if n_lower > SEQ_LEN // 2:
            n_repeat += 1
            continue
        kept.append(w.upper())

    print(
        f"kept {len(kept):,} | oob={n_oob:,} bad_char={n_bad_char:,} repeat={n_repeat:,}",
        file=sys.stderr,
    )

    if len(kept) < N_SEQS:
        raise SystemExit(f"Only {len(kept):,} survivors; need {N_SEQS:,}.")

    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(kept), size=N_SEQS, replace=False)
    sample = [kept[i] for i in idx]

    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in sample[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(sample))
        f.write("\n")

    print(f"Wrote {len(sample):,} sequences to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
