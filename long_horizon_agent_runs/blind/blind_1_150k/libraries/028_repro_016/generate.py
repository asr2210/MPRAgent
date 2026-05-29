"""Experiment 016 — Multi-source: cCRE (mixed strand) + FANTOM5 + motif.

Tests whether FANTOM5 CAGE-defined enhancers (alternative annotation
methodology) add information beyond ENCODE cCREs when included as part
of the diversity component.

Composition:
- 67,500 cCREs (forward only), seed=15
- 67,500 DIFFERENT cCREs in reverse complement, seed=15
- 7,500 FANTOM5 enhancers, 200bp centered, ACGT + softmask-filtered, seed=15
- 7,500 motif-embedded synthetic, seed=15
- Total: 150k.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

N_FWD = 67_500
N_RC = 67_500
N_FANTOM = 7_500
N_SYNTH = 7_500
N_SEQS = 150_000
N_MOTIFS_PER_SEQ = 2
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 27
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
FANTOM_PATH = os.path.join(DATA_DIR, "fantom5_hg38_enhancers.bed")
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


def passes(w: str) -> bool:
    valid = set("ACGTacgt")
    if any(c not in valid for c in w):
        return False
    n_lower = sum(1 for c in w if c.islower())
    return n_lower <= SEQ_LEN // 2


def build_bed_window_pool(
    bed_path: str, chrom_seq: dict[str, str]
) -> list[str]:
    kept: list[str] = []
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            mid = (start + end) // 2
            seq = chrom_seq.get(chrom)
            if seq is None:
                continue
            lo, hi = mid - HALF, mid + HALF
            if lo < 0 or hi > len(seq):
                continue
            w = seq[lo:hi]
            if not passes(w):
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

    # Load all needed chromosomes (cCREs cover all; FANTOM5 mostly the same).
    chroms_needed: set[str] = set()
    for path in (BED_PATH, FANTOM_PATH):
        with open(path) as f:
            for line in f:
                chroms_needed.add(line.split("\t", 1)[0])
    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        p = os.path.join(DATA_DIR, f"{c}.fa")
        if os.path.exists(p):
            chrom_seq[c] = read_fasta_seq(p)

    print("Building cCRE pool...", file=sys.stderr)
    ccre_pool = build_bed_window_pool(BED_PATH, chrom_seq)
    print(f"  cCRE pool: {len(ccre_pool):,}", file=sys.stderr)

    print("Building FANTOM5 enhancer pool...", file=sys.stderr)
    fantom_pool = build_bed_window_pool(FANTOM_PATH, chrom_seq)
    print(f"  FANTOM5 pool: {len(fantom_pool):,}", file=sys.stderr)
    assert len(fantom_pool) >= N_FANTOM

    n_needed_ccre = N_FWD + N_RC
    assert len(ccre_pool) >= n_needed_ccre

    ccre_idx = rng.choice(len(ccre_pool), size=n_needed_ccre, replace=False)
    fwd = [ccre_pool[i] for i in ccre_idx[:N_FWD]]
    rc = [revcomp(ccre_pool[i]) for i in ccre_idx[N_FWD:]]

    fantom_idx = rng.choice(len(fantom_pool), size=N_FANTOM, replace=False)
    fantom = [fantom_pool[i] for i in fantom_idx]

    motifs = load_jaspar_consensus(JASPAR_PATH)
    synth = make_motif_embedded(N_SYNTH, motifs, rng)

    combined = fwd + rc + fantom + synth
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote {N_SEQS:,} ({N_FWD} fwd + {N_RC} RC + {N_FANTOM} fantom + {N_SYNTH} motif)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
