"""Experiment 008 — 90/10 cCRE + non-cCRE genomic tiles.

Compares biological diversity (genomic windows outside any cCRE) against
synthetic motif-embedded random as the 10% diversity source in the winning
90/10 hybrid. Direct comparison to experiment 006 (mean 0.888).

Composition:
- 135,000 cCRE-centered 200bp windows (same pool & filter), seed=7
- 15,000 non-cCRE genomic 200bp tiles: hg38 chr1, 7, 14, 19-22 tiled
  non-overlap, filter (ACGT + <50% softmasked), exclude any window whose
  centre is within 500 bp of an ENCODE cCRE midpoint (so the diversity
  component is genuinely non-cCRE), sampled uniformly with seed=7
- Concatenated, shuffled, written.
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
N_CCRE = 135_000
N_NONCCRE = 15_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
CCRE_EXCLUSION_DISTANCE = 500  # bp from nearest cCRE midpoint to keep
SEED = 7
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
NONCCRE_CHROMS = ("chr1", "chr7", "chr14", "chr19", "chr20", "chr21", "chr22")
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


def build_ccre_pool(
    chrom_seq: dict[str, str], ccre_mids_by_chrom: dict[str, list[int]]
) -> list[str]:
    valid = set("ACGTacgt")
    kept: list[str] = []
    for chrom, mids in ccre_mids_by_chrom.items():
        seq = chrom_seq.get(chrom)
        if seq is None:
            continue
        for mid in mids:
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


def build_nonccre_tile_pool(
    chrom_seq: dict[str, str],
    ccre_mids_by_chrom: dict[str, list[int]],
) -> list[str]:
    """Tile target chromosomes; drop tiles near any cCRE midpoint."""
    valid = set("ACGTacgt")
    kept: list[str] = []
    for c in NONCCRE_CHROMS:
        seq = chrom_seq.get(c)
        if seq is None:
            continue
        # Sort cCRE mids for binary search.
        mids = sorted(ccre_mids_by_chrom.get(c, []))
        mids_arr = np.array(mids, dtype=np.int64)
        n_tiles = len(seq) // SEQ_LEN
        for i in range(n_tiles):
            lo = i * SEQ_LEN
            hi = lo + SEQ_LEN
            w = seq[lo:hi]
            if any(ch not in valid for ch in w):
                continue
            if sum(1 for ch in w if ch.islower()) > SEQ_LEN // 2:
                continue
            mid = (lo + hi) // 2
            # nearest cCRE midpoint distance check
            if len(mids_arr) > 0:
                idx = int(np.searchsorted(mids_arr, mid))
                near = []
                if idx < len(mids_arr):
                    near.append(abs(int(mids_arr[idx]) - mid))
                if idx > 0:
                    near.append(abs(int(mids_arr[idx - 1]) - mid))
                if near and min(near) <= CCRE_EXCLUSION_DISTANCE:
                    continue
            kept.append(w.upper())
    return kept


def main() -> None:
    rng = np.random.default_rng(SEED)

    print("Loading cCRE bed and chromosomes...", file=sys.stderr)
    ccre_mids_by_chrom: dict[str, list[int]] = {}
    chroms_needed: set[str] = set()
    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            mid = (int(parts[1]) + int(parts[2])) // 2
            ccre_mids_by_chrom.setdefault(chrom, []).append(mid)
            chroms_needed.add(chrom)
    chroms_needed.update(NONCCRE_CHROMS)

    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        path = os.path.join(DATA_DIR, f"{c}.fa")
        if os.path.exists(path):
            chrom_seq[c] = read_fasta_seq(path)

    print("Building cCRE pool...", file=sys.stderr)
    ccre_pool = build_ccre_pool(chrom_seq, ccre_mids_by_chrom)
    print(f"  cCRE pool: {len(ccre_pool):,}", file=sys.stderr)
    assert len(ccre_pool) >= N_CCRE
    ccre_idx = rng.choice(len(ccre_pool), size=N_CCRE, replace=False)
    ccre_sample = [ccre_pool[i] for i in ccre_idx]

    print("Building non-cCRE genomic tile pool...", file=sys.stderr)
    nonccre_pool = build_nonccre_tile_pool(chrom_seq, ccre_mids_by_chrom)
    print(f"  non-cCRE tile pool: {len(nonccre_pool):,}", file=sys.stderr)
    assert len(nonccre_pool) >= N_NONCCRE
    nonccre_idx = rng.choice(len(nonccre_pool), size=N_NONCCRE, replace=False)
    nonccre_sample = [nonccre_pool[i] for i in nonccre_idx]

    combined = ccre_sample + nonccre_sample
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote {N_SEQS:,} ({N_CCRE:,} cCRE + {N_NONCCRE:,} non-cCRE genome) to {OUT_PATH}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
