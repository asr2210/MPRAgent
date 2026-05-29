"""Experiment 009 — 90/10 cCRE + dense motif-embedded (5 motifs/seq).

Direct refinement of exp 006: same 135k cCRE + 15k synthetic split, but
with 5 JASPAR consensus motifs per synthetic 200bp sequence instead of 2.
Tests whether denser motif content within the diversity component improves
mean_r without losing the eval_08 random-like bonus.

If positions would overlap, retries up to 50 times; if still no free slot,
just stamps over.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

N_SEQS = 150_000
N_CCRE = 135_000
N_SYNTH = 15_000
N_MOTIFS_PER_SEQ = 5
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 8
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
JASPAR_PATH = os.path.join(DATA_DIR, "jaspar2024_core_vert_pfms.txt")
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


_BASE_ROW_RE = re.compile(r"^([ACGT])\s*\[([^\]]+)\]\s*$")


def load_jaspar_consensus(path: str) -> list[str]:
    consensus_list: list[str] = []
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue
        rows: list[list[float]] = [[], [], [], []]
        for j in range(4):
            m = _BASE_ROW_RE.match(lines[i + 1 + j])
            if not m:
                raise ValueError(f"bad JASPAR row at line {i+1+j}")
            base = m.group(1)
            counts = [float(x) for x in m.group(2).split()]
            rows[base_to_idx[base]] = counts
        L = len(rows[0])
        if not all(len(r) == L for r in rows):
            raise ValueError(f"non-rectangular PFM at line {i}")
        counts = np.array(rows)
        argmax = np.argmax(counts, axis=0)
        consensus = "".join("ACGT"[k] for k in argmax)
        if 4 <= len(consensus) <= SEQ_LEN:
            consensus_list.append(consensus)
        i += 5
    return consensus_list


def make_dense_motif_embedded(
    n: int, motifs: list[str], rng: np.random.Generator
) -> list[str]:
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    n_motifs = len(motifs)
    motif_choice = rng.integers(0, n_motifs, size=(n, N_MOTIFS_PER_SEQ))
    out: list[str] = []
    for i in range(n):
        bg = chars[i].copy()
        occupied = np.zeros(SEQ_LEN, dtype=bool)
        for k in range(N_MOTIFS_PER_SEQ):
            m = motifs[int(motif_choice[i, k])]
            L = len(m)
            max_start = SEQ_LEN - L
            # try up to 50 times to find a non-overlapping placement
            placed = False
            for _ in range(50):
                pos = int(rng.integers(0, max_start + 1))
                if not occupied[pos : pos + L].any():
                    bg[pos : pos + L] = np.array(list(m))
                    occupied[pos : pos + L] = True
                    placed = True
                    break
            if not placed:
                # fallback: stamp anyway (overwrites a previous motif partially)
                pos = int(rng.integers(0, max_start + 1))
                bg[pos : pos + L] = np.array(list(m))
        out.append("".join(bg.tolist()))
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)

    ccre_pool = build_ccre_pool()
    print(f"cCRE pool: {len(ccre_pool):,}", file=sys.stderr)
    ccre_idx = rng.choice(len(ccre_pool), size=N_CCRE, replace=False)
    ccre_sample = [ccre_pool[i] for i in ccre_idx]

    motifs = load_jaspar_consensus(JASPAR_PATH)
    print(f"motifs: {len(motifs)}", file=sys.stderr)
    synth = make_dense_motif_embedded(N_SYNTH, motifs, rng)
    print(f"synth: {len(synth):,}", file=sys.stderr)

    combined = ccre_sample + synth
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote {N_SEQS:,} ({N_CCRE:,} cCRE + {N_SYNTH:,} dense motif-embedded) to {OUT_PATH}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
