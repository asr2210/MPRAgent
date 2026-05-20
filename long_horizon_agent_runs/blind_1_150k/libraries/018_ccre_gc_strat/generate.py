"""Experiment 018 — GC-stratified cCRE selection.

Forces a uniform GC distribution over the cCRE pool, then applies the
winning 015 recipe (mixed-strand + 15k motif).

Composition:
- 67,500 cCREs (forward) drawn uniformly across 30 GC bins, seed=17
- 67,500 different cCREs (reverse complement), drawn the same way
- 15,000 motif-embedded synthetic
- Total 150k.

Tests whether stratifying cCRE selection on a continuous axis (GC) helps
generalisation (vs class-stratification which hurt in exp 007).
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

N_FWD = 67_500
N_RC = 67_500
N_SYNTH = 15_000
N_SEQS = 150_000
N_BINS = 30
N_MOTIFS_PER_SEQ = 2
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 17
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


def gc_content(s: str) -> float:
    return sum(1 for c in s if c in "GC") / len(s)


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


def main() -> None:
    rng = np.random.default_rng(SEED)
    ccre_pool = build_ccre_pool()
    print(f"cCRE pool: {len(ccre_pool):,}", file=sys.stderr)

    gcs = np.array([gc_content(s) for s in ccre_pool])
    bin_edges = np.linspace(gcs.min() - 1e-9, gcs.max() + 1e-9, N_BINS + 1)
    bin_idx = np.digitize(gcs, bin_edges) - 1  # 0..N_BINS-1
    print("GC bins:", file=sys.stderr)
    for b in range(N_BINS):
        n_in_bin = int((bin_idx == b).sum())
        lo, hi = bin_edges[b], bin_edges[b + 1]
        print(f"  bin {b:2d}: GC [{lo:.3f}, {hi:.3f}) — {n_in_bin:,}", file=sys.stderr)

    per_bin_fwd = N_FWD // N_BINS  # 2250
    per_bin_rc = N_RC // N_BINS    # 2250
    # Distribute the remainder over the highest-count bins:
    fwd_remainder = N_FWD - per_bin_fwd * N_BINS
    rc_remainder = N_RC - per_bin_rc * N_BINS

    # Per-bin index lists
    bin_to_idx: dict[int, np.ndarray] = {b: np.where(bin_idx == b)[0] for b in range(N_BINS)}

    fwd_selected_global: list[int] = []
    rc_selected_global: list[int] = []
    for b in range(N_BINS):
        idxs = bin_to_idx[b]
        need_fwd = per_bin_fwd + (1 if b < fwd_remainder else 0)
        need_rc = per_bin_rc + (1 if b < rc_remainder else 0)
        if len(idxs) >= need_fwd + need_rc:
            picks = rng.choice(idxs, size=need_fwd + need_rc, replace=False)
        else:
            picks = rng.choice(idxs, size=need_fwd + need_rc, replace=True)
        fwd_selected_global.extend(picks[:need_fwd].tolist())
        rc_selected_global.extend(picks[need_fwd:].tolist())

    fwd = [ccre_pool[i] for i in fwd_selected_global]
    rc = [revcomp(ccre_pool[i]) for i in rc_selected_global]
    print(f"selected {len(fwd)} fwd + {len(rc)} RC", file=sys.stderr)

    motifs = load_jaspar_consensus(JASPAR_PATH)
    synth = make_motif_embedded(N_SYNTH, motifs, rng)

    combined = fwd + rc + synth
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(f"Wrote {N_SEQS:,} sequences to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
