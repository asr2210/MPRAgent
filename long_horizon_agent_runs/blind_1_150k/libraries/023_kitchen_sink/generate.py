"""Experiment 023 — Kitchen-sink: combine all winning levers.

Tests whether stacking the winning ingredients gives a compound benefit
above the ±0.003 training-noise floor.

Composition:
- 60,000 cCRE-fwd + 60,000 different cCRE-RC (120k unique cCREs, mixed strand)
- 7,500 FANTOM5-fwd + 7,500 different FANTOM5-RC (15k FANTOM5, mixed strand)
- 10,000 motif-embedded + 5,000 uniform random synthetic (15k synth, 2 sources)
- Total 150k, seed=22
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

N_CCRE_FWD = 60_000
N_CCRE_RC = 60_000
N_FANTOM_FWD = 7_500
N_FANTOM_RC = 7_500
N_SYNTH_MOTIF = 10_000
N_SYNTH_RAND = 5_000
N_SEQS = 150_000
N_MOTIFS_PER_SEQ = 2
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 22
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
    return sum(1 for c in w if c.islower()) <= SEQ_LEN // 2


def build_bed_pool(bed_path: str, chrom_seq: dict[str, str]) -> list[str]:
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


def make_uniform_random(n: int, rng: np.random.Generator) -> list[str]:
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]
    return chars.view(f"<U{SEQ_LEN}").ravel().tolist()


def main() -> None:
    rng = np.random.default_rng(SEED)

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

    print("cCRE pool...", file=sys.stderr)
    ccre_pool = build_bed_pool(BED_PATH, chrom_seq)
    print(f"  {len(ccre_pool):,}", file=sys.stderr)
    print("FANTOM5 pool...", file=sys.stderr)
    fantom_pool = build_bed_pool(FANTOM_PATH, chrom_seq)
    print(f"  {len(fantom_pool):,}", file=sys.stderr)

    n_ccre_needed = N_CCRE_FWD + N_CCRE_RC
    n_fantom_needed = N_FANTOM_FWD + N_FANTOM_RC
    ccre_idx = rng.choice(len(ccre_pool), size=n_ccre_needed, replace=False)
    fantom_idx = rng.choice(len(fantom_pool), size=n_fantom_needed, replace=False)

    ccre_fwd = [ccre_pool[i] for i in ccre_idx[:N_CCRE_FWD]]
    ccre_rc = [revcomp(ccre_pool[i]) for i in ccre_idx[N_CCRE_FWD:]]
    fantom_fwd = [fantom_pool[i] for i in fantom_idx[:N_FANTOM_FWD]]
    fantom_rc = [revcomp(fantom_pool[i]) for i in fantom_idx[N_FANTOM_FWD:]]

    motifs = load_jaspar_consensus(JASPAR_PATH)
    synth_motif = make_motif_embedded(N_SYNTH_MOTIF, motifs, rng)
    synth_rand = make_uniform_random(N_SYNTH_RAND, rng)

    combined = ccre_fwd + ccre_rc + fantom_fwd + fantom_rc + synth_motif + synth_rand
    assert len(combined) == N_SEQS
    rng.shuffle(combined)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")
    print(
        f"Wrote {N_SEQS:,} (cCRE {N_CCRE_FWD + N_CCRE_RC}, FANTOM {n_fantom_needed}, "
        f"motif {N_SYNTH_MOTIF}, random {N_SYNTH_RAND})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
