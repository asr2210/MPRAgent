"""Experiment 004 — 50/50 hybrid: cCRE-centered + uniform random.

Tests the diversity hypothesis suggested by eval_08 in experiment 002 and
the partial recovery in experiment 003. Random had a unique strength on
eval_08 that cCREs only partially closed; this experiment asks whether
combining the two sources strictly beats either alone.

Composition:
- 75,000 cCRE-centered 200bp windows (same pool & filter as exp 003), seed=3
- 75,000 uniform random 200bp ACGT sequences, seed=3
- Concatenated and shuffled (seed=3) before writing.

Comparison:
- Experiment 001 (random):  mean_r 0.820  (eval_08 = 0.908)
- Experiment 002 (genome):  mean_r 0.873  (eval_08 = 0.868)
- Experiment 003 (cCRE):    mean_r 0.886  (eval_08 = 0.892)
- Experiment 004 (hybrid):  ?
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
N_CCRE = 75_000
N_RAND = 75_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 3
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
ALPHABET = np.array(list("ACGT"))
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


def build_ccre_pool() -> list[str]:
    chroms_needed: set[str] = set()
    cres: list[tuple[str, int]] = []
    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            mid = (int(parts[1]) + int(parts[2])) // 2
            cres.append((chrom, mid))
            chroms_needed.add(chrom)

    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        path = os.path.join(DATA_DIR, f"{c}.fa")
        if os.path.exists(path):
            chrom_seq[c] = read_fasta_seq(path)

    valid = set("ACGTacgt")
    kept: list[str] = []
    for chrom, mid in cres:
        seq = chrom_seq.get(chrom)
        if seq is None:
            continue
        lo, hi = mid - HALF, mid + HALF
        if lo < 0 or hi > len(seq):
            continue
        w = seq[lo:hi]
        if any(c not in valid for c in w):
            continue
        if sum(1 for c in w if c.islower()) > SEQ_LEN // 2:
            continue
        kept.append(w.upper())
    return kept


def sample_random(n: int, rng: np.random.Generator) -> list[str]:
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    return chars.view(f"<U{SEQ_LEN}").ravel().tolist()


def main() -> None:
    rng = np.random.default_rng(SEED)

    print("Building cCRE pool...", file=sys.stderr)
    ccre_pool = build_ccre_pool()
    print(f"  cCRE pool: {len(ccre_pool):,} survivors", file=sys.stderr)
    assert len(ccre_pool) >= N_CCRE

    ccre_idx = rng.choice(len(ccre_pool), size=N_CCRE, replace=False)
    ccre_sample = [ccre_pool[i] for i in ccre_idx]

    rand_sample = sample_random(N_RAND, rng)

    combined = ccre_sample + rand_sample
    assert len(combined) == N_SEQS
    rng.shuffle(combined)

    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")

    print(f"Wrote {N_SEQS:,} sequences ({N_CCRE:,} cCRE + {N_RAND:,} random) to {OUT_PATH}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
