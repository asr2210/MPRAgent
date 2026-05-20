"""Experiment 011 — cCRE + mutated cCRE pairs.

For each of 75,000 unique cCRE-centered windows, generate a "mutated"
counterpart by randomly substituting 10% of bases (20 of 200) to a
different ACGT base. Library = 75k originals + 75k mutated = 150k.

Tests whether data augmentation by small perturbations gives the model
contrastive training signal — paired (sequence, activity) and (mutated,
activity') examples — that lifts performance beyond what 150k unique
cCRE sequences alone provide (exp 003 mean=0.886).
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_PAIRS = 75_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
N_MUTATE = SEQ_LEN // 10  # 20 bases per 200
SEED = 10
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")

ALPHABET = np.array(list("ACGT"))
BASE_TO_IDX = {b: i for i, b in enumerate("ACGT")}


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


def mutate(seq: str, n: int, rng: np.random.Generator) -> str:
    """Substitute n random positions with a different ACGT base."""
    arr = np.array(list(seq))
    positions = rng.choice(SEQ_LEN, size=n, replace=False)
    for p in positions:
        orig = arr[p]
        idx = BASE_TO_IDX[orig]
        # pick a different base
        new_idx = int(rng.integers(0, 3))
        if new_idx >= idx:
            new_idx += 1
        arr[p] = "ACGT"[new_idx]
    return "".join(arr.tolist())


def main() -> None:
    rng = np.random.default_rng(SEED)
    ccre_pool = build_ccre_pool()
    print(f"cCRE pool: {len(ccre_pool):,}", file=sys.stderr)
    assert len(ccre_pool) >= N_PAIRS

    idx = rng.choice(len(ccre_pool), size=N_PAIRS, replace=False)
    originals = [ccre_pool[i] for i in idx]
    mutated = [mutate(s, N_MUTATE, rng) for s in originals]

    # Sanity: original and mutated differ by exactly N_MUTATE positions.
    diffs = sum(a != b for a, b in zip(originals[0], mutated[0]))
    print(f"sanity: first pair Hamming distance = {diffs} (expected {N_MUTATE})", file=sys.stderr)
    assert diffs == N_MUTATE

    combined = originals + mutated
    assert len(combined) == 2 * N_PAIRS == 150_000
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote 150,000 sequences ({N_PAIRS} unique cCREs + their 10%-mutated copies) "
        f"to {OUT_PATH}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
