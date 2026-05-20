"""Experiment 005 — 50/50 hybrid: cCRE-centered + motif-embedded synthetic.

Replaces the uniform-random half of experiment 004 with a more focused
diversity source: random 200bp backgrounds with two JASPAR motif consensus
sequences embedded at random positions. Tests whether 'motif content' is
the active ingredient of random's contribution to the hybrid (eval_08 win).

Composition:
- 75,000 cCRE-centered 200bp windows, seed=4 (same pool as exp 003/004).
- 75,000 synthetic 200bp = random ACGT background + 2 JASPAR motifs embedded
  (consensus base per PFM position) at random positions, seed=4. Motifs are
  sampled uniformly from the JASPAR 2024 vertebrate non-redundant set
  (~880 motifs).
- Concatenated, shuffled, written.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

N_SEQS = 150_000
N_CCRE = 75_000
N_SYNTH = 75_000
N_MOTIFS_PER_SEQ = 2
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 4
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
    """Parse JASPAR PFM file and return list of consensus DNA strings.

    Each motif occupies a header line ('>') followed by 4 base-count rows
    in order A, C, G, T. Consensus = argmax base at each position.
    """
    consensus_list: list[str] = []
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue
        # next 4 lines are the matrix rows
        rows: list[list[int]] = [[], [], [], []]
        for j in range(4):
            m = _BASE_ROW_RE.match(lines[i + 1 + j])
            if not m:
                raise ValueError(f"bad JASPAR row at line {i+1+j}: {lines[i+1+j]!r}")
            base = m.group(1)
            counts = [float(x) for x in m.group(2).split()]
            rows[base_to_idx[base]] = counts
        # All four rows must be same length (motif length).
        L = len(rows[0])
        if not all(len(r) == L for r in rows):
            raise ValueError(f"non-rectangular PFM at line {i}")
        # argmax per column → consensus.
        counts = np.array(rows)  # (4, L)
        argmax = np.argmax(counts, axis=0)  # (L,)
        consensus = "".join("ACGT"[k] for k in argmax)
        if 4 <= len(consensus) <= SEQ_LEN:  # sanity limit
            consensus_list.append(consensus)
        i += 5
    return consensus_list


def make_motif_embedded(n: int, motifs: list[str], rng: np.random.Generator) -> list[str]:
    # Random ACGT backgrounds.
    idx = rng.integers(0, 4, size=(n, SEQ_LEN), dtype=np.uint8)
    chars = ALPHABET[idx]  # (n, SEQ_LEN)
    # Embed N_MOTIFS_PER_SEQ motifs per sequence at random positions.
    n_motifs = len(motifs)
    motif_choice = rng.integers(0, n_motifs, size=(n, N_MOTIFS_PER_SEQ))
    out: list[str] = []
    for i in range(n):
        bg = chars[i].copy()
        for k in range(N_MOTIFS_PER_SEQ):
            m = motifs[int(motif_choice[i, k])]
            L = len(m)
            max_start = SEQ_LEN - L
            pos = int(rng.integers(0, max_start + 1))
            bg[pos : pos + L] = np.array(list(m))
        out.append("".join(bg.tolist()))
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)

    print("Building cCRE pool...", file=sys.stderr)
    ccre_pool = build_ccre_pool()
    print(f"  cCRE pool: {len(ccre_pool):,} survivors", file=sys.stderr)
    assert len(ccre_pool) >= N_CCRE
    ccre_idx = rng.choice(len(ccre_pool), size=N_CCRE, replace=False)
    ccre_sample = [ccre_pool[i] for i in ccre_idx]

    print("Loading JASPAR motifs...", file=sys.stderr)
    motifs = load_jaspar_consensus(JASPAR_PATH)
    print(f"  loaded {len(motifs)} motif consensus sequences", file=sys.stderr)
    lengths = [len(m) for m in motifs]
    print(
        f"  motif length: min={min(lengths)} max={max(lengths)} "
        f"mean={np.mean(lengths):.1f}",
        file=sys.stderr,
    )

    print("Building synthetic motif-embedded backgrounds...", file=sys.stderr)
    synth = make_motif_embedded(N_SYNTH, motifs, rng)
    print(f"  built {len(synth):,} synthetic sequences", file=sys.stderr)

    combined = ccre_sample + synth
    assert len(combined) == N_SEQS
    rng.shuffle(combined)

    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in combined[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(combined))
        f.write("\n")

    print(
        f"Wrote {N_SEQS:,} sequences ({N_CCRE:,} cCRE + {N_SYNTH:,} motif-embedded) to {OUT_PATH}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
