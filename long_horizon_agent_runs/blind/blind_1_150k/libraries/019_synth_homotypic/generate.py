"""Experiment 019 — Synthetic with homotypic motif clusters.

Tests whether biologically-motivated synthetic architecture (3 copies of
the same JASPAR motif spaced ~30bp apart, mimicking homotypic TF clusters)
beats random-placement motif synthetic as the diversity component.

Composition:
- 67,500 cCRE-fwd + 67,500 cCRE-RC of different (135k unique, mixed strand)
- 15,000 homotypic-cluster synthetic: random 200bp ACGT background with 3
  copies of the SAME JASPAR motif at positions ~50, 100, 150
- Total 150k.
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
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 18
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


def make_homotypic_cluster_synth(
    n: int, motifs: list[str], rng: np.random.Generator
) -> list[str]:
    """Each sequence has 3 copies of one motif, evenly spaced."""
    # Only use motifs of length <= 30bp (need 3 to fit at positions ~50, 100, 150
    # with 30bp gaps).
    short_motifs = [m for m in motifs if len(m) <= 30]
    print(f"  using {len(short_motifs)} motifs ≤ 30bp for homotypic clusters", file=sys.stderr)

    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    motif_choice = rng.integers(0, len(short_motifs), size=n)
    # base positions: 50, 100, 150 (centres for 30bp-spaced clusters)
    # but adjust to fit motif length: position = anchor - L/2
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
    synth = make_homotypic_cluster_synth(N_SYNTH, motifs, rng)

    combined = fwd + rc + synth
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote {N_SEQS:,} ({N_FWD} fwd + {N_RC} RC + {N_SYNTH} homotypic-cluster synth)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
