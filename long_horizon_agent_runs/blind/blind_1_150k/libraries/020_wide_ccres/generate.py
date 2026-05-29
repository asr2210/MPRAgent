"""Experiment 020 — Wide cCREs (>= 300bp) only.

Filters cCREs to width >= 300bp under the hypothesis that "wider" cCREs
are more complex regulatory regions containing richer architecture per
sequence. Combined with 015 recipe (mixed-strand + 15k motif).

Pool: 488k cCREs ≥300bp out of 1.06M total (~46%). After standard filter:
expect ~320k.

Composition:
- 67,500 wide cCREs (forward), seed=19
- 67,500 DIFFERENT wide cCREs (reverse complement), seed=19
- 15,000 motif-embedded synthetic, seed=19
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
N_MOTIFS_PER_SEQ = 2
MIN_CCRE_WIDTH = 300
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 19
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


def build_wide_pool() -> list[str]:
    chroms_needed: set[str] = set()
    cres: list[tuple[str, int]] = []
    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if end - start < MIN_CCRE_WIDTH:
                continue
            mid = (start + end) // 2
            cres.append((chrom, mid))
            chroms_needed.add(chrom)
    print(f"cCREs with width >= {MIN_CCRE_WIDTH}bp: {len(cres):,}", file=sys.stderr)
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


def main() -> None:
    rng = np.random.default_rng(SEED)
    wide_pool = build_wide_pool()
    print(f"wide cCRE pool (after filter): {len(wide_pool):,}", file=sys.stderr)
    n_needed = N_FWD + N_RC
    assert len(wide_pool) >= n_needed
    selection = rng.choice(len(wide_pool), size=n_needed, replace=False)
    fwd = [wide_pool[i] for i in selection[:N_FWD]]
    rc = [revcomp(wide_pool[i]) for i in selection[N_FWD:]]

    motifs = load_jaspar_consensus(JASPAR_PATH)
    synth = make_motif_embedded(N_SYNTH, motifs, rng)

    combined = fwd + rc + synth
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(f"Wrote {N_SEQS:,} sequences ({N_FWD} fwd + {N_RC} RC + {N_SYNTH} motif)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
