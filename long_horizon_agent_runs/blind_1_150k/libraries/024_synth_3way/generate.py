"""Experiment 024 — 3-source synthetic split.

Same 015 cCRE bulk but with synthetic split across three sources:
motif-embedded (random placement), uniform random, and homotypic clusters.
Tests if synthetic source diversity has any compounding effect.

Composition:
- 67,500 cCRE-fwd + 67,500 different cCRE-RC (135k unique, mixed strand)
- 5,000 motif-embedded (2 motifs random placement)
- 5,000 uniform random ACGT
- 5,000 homotypic clusters (3x same motif at ~50/100/150)
- Total 150k, seed=23
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

N_FWD = 67_500
N_RC = 67_500
N_MOTIF = 5_000
N_RAND = 5_000
N_CLUSTER = 5_000
N_SEQS = 150_000
N_MOTIFS_PER_SEQ = 2
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 23
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
JASPAR_PATH = os.path.join(DATA_DIR, "jaspar2024_core_vert_pfms.txt")
ALPHABET = np.array(list("ACGT"))
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")

RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(RC_TABLE)[::-1]


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
    out: list[str] = []
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
                raise ValueError(f"bad row {i+1+j}")
            base = m.group(1)
            counts = [float(x) for x in m.group(2).split()]
            rows[base_to_idx[base]] = counts
        L = len(rows[0])
        if not all(len(r) == L for r in rows):
            raise ValueError("non-rectangular PFM")
        counts = np.array(rows)
        consensus = "".join("ACGT"[k] for k in np.argmax(counts, axis=0))
        if 4 <= len(consensus) <= SEQ_LEN:
            out.append(consensus)
        i += 5
    return out


def make_motif_embedded(n: int, motifs: list[str], rng: np.random.Generator) -> list[str]:
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    nm = len(motifs)
    motif_choice = rng.integers(0, nm, size=(n, N_MOTIFS_PER_SEQ))
    out: list[str] = []
    for i in range(n):
        bg = chars[i].copy()
        for k in range(N_MOTIFS_PER_SEQ):
            m = motifs[int(motif_choice[i, k])]
            L = len(m)
            pos = int(rng.integers(0, SEQ_LEN - L + 1))
            bg[pos : pos + L] = np.array(list(m))
        out.append("".join(bg.tolist()))
    return out


def make_uniform_random(n: int, rng: np.random.Generator) -> list[str]:
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    return chars.view(f"<U{SEQ_LEN}").ravel().tolist()


def make_homotypic_clusters(n: int, motifs: list[str], rng: np.random.Generator) -> list[str]:
    short_motifs = [m for m in motifs if len(m) <= 30]
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    motif_choice = rng.integers(0, len(short_motifs), size=n)
    anchors = [50, 100, 150]
    out: list[str] = []
    for i in range(n):
        bg = chars[i].copy()
        m = short_motifs[int(motif_choice[i])]
        L = len(m)
        for anchor in anchors:
            pos = max(0, min(SEQ_LEN - L, anchor - L // 2))
            bg[pos : pos + L] = np.array(list(m))
        out.append("".join(bg.tolist()))
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    ccre_pool = build_ccre_pool()
    print(f"cCRE pool: {len(ccre_pool):,}", file=sys.stderr)
    n_needed = N_FWD + N_RC
    selection = rng.choice(len(ccre_pool), size=n_needed, replace=False)
    fwd = [ccre_pool[i] for i in selection[:N_FWD]]
    rc = [revcomp(ccre_pool[i]) for i in selection[N_FWD:]]

    motifs = load_jaspar_consensus(JASPAR_PATH)
    motif_seqs = make_motif_embedded(N_MOTIF, motifs, rng)
    rand_seqs = make_uniform_random(N_RAND, rng)
    cluster_seqs = make_homotypic_clusters(N_CLUSTER, motifs, rng)

    combined = fwd + rc + motif_seqs + rand_seqs + cluster_seqs
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote {N_SEQS:,} ({N_FWD} fwd + {N_RC} RC + {N_MOTIF} motif + {N_RAND} random "
        f"+ {N_CLUSTER} cluster)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
