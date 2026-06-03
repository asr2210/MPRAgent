"""
001_motif_cocktail
==================
50,000 200bp synthetic sequences. Each = random ~50% GC background with
4-8 TF motif instances inserted at random non-overlapping positions.
Motifs sampled from JASPAR 2024 CORE PFMs (2346 motifs).

Hypothesis: planted motif content alone is enough to beat the pure-random
synthetic baseline (synth_oracle eval_01 = 0.6840). Whether it beats
DHS-derived libraries (dhs_topic = 0.7232) tests whether genomic context
matters above and beyond motif density.

Output: sequences_0.txt (50000 lines, 200 chars from {A,C,G,T})
"""
import os
import random
import sys
from pathlib import Path

import numpy as np

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
MIN_MOTIFS = 4
MAX_MOTIFS = 8
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
JASPAR_PATH = DATA_DIR / "jaspar2024_core.txt"


def parse_jaspar_pfms(path: Path):
    """Parse JASPAR PFM text format. Returns list of (name, pfm) where pfm is
    a (4, L) numpy array of column-normalized probabilities (A,C,G,T)."""
    pfms = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith(">"):
            name = line[1:].strip()
            rows = []
            for j in range(1, 5):
                row = lines[i + j].strip()
                # format: "A  [   4   19   ...   ]"
                start = row.index("[") + 1
                end = row.index("]")
                vals = [float(x) for x in row[start:end].split()]
                rows.append(vals)
            mat = np.array(rows, dtype=np.float64)  # 4 x L
            # column-normalize with small pseudocount
            col_sums = mat.sum(axis=0, keepdims=True)
            col_sums[col_sums == 0] = 1.0
            mat = (mat + 0.01) / (col_sums + 0.04)
            pfms.append((name, mat))
            i += 5
        else:
            i += 1
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    """Sample one instance of the motif by drawing each position from its PFM."""
    bases = np.array(["A", "C", "G", "T"])
    L = pfm.shape[1]
    out = []
    for c in range(L):
        out.append(rng.choice(bases, p=pfm[:, c]))
    return "".join(out)


def random_bg(rng: np.random.Generator, length: int) -> list:
    """Random 50% GC background as list of chars (so we can splice motifs in)."""
    bases = np.array(["A", "C", "G", "T"])
    return list(rng.choice(bases, size=length))


def plant_motifs(seq: list, motifs: list, rng: np.random.Generator) -> str:
    """Plant each motif string at a random non-overlapping position."""
    L = len(seq)
    occupied = []  # list of (start, end) intervals
    for m in motifs:
        ml = len(m)
        if ml >= L:
            continue
        # try up to 50 positions to find a non-overlap
        for _ in range(50):
            start = int(rng.integers(0, L - ml + 1))
            end = start + ml
            if all(end <= a or start >= b for a, b in occupied):
                break
        else:
            # Could not find an empty slot; overwrite anyway
            start = int(rng.integers(0, L - ml + 1))
            end = start + ml
        occupied.append((start, end))
        for k, ch in enumerate(m):
            seq[start + k] = ch
    return "".join(seq)


def main():
    rng = np.random.default_rng(SEED)
    random.seed(SEED)

    print(f"Parsing JASPAR PFMs from {JASPAR_PATH}", file=sys.stderr)
    pfms = parse_jaspar_pfms(JASPAR_PATH)
    print(f"Parsed {len(pfms)} motifs", file=sys.stderr)

    # Filter to reasonable length motifs
    pfms = [(n, m) for n, m in pfms if 6 <= m.shape[1] <= 20]
    print(f"After length filter (6<=L<=20): {len(pfms)} motifs", file=sys.stderr)

    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for i in range(N_SEQS):
            n_motifs = int(rng.integers(MIN_MOTIFS, MAX_MOTIFS + 1))
            chosen_idx = rng.choice(len(pfms), size=n_motifs, replace=True)
            motif_instances = [sample_from_pfm(pfms[k][1], rng) for k in chosen_idx]
            bg = random_bg(rng, SEQ_LEN)
            seq = plant_motifs(bg, motif_instances, rng)
            assert len(seq) == SEQ_LEN
            assert set(seq) <= {"A", "C", "G", "T"}
            f.write(seq + "\n")
            if (i + 1) % 10_000 == 0:
                print(f"  {i+1}/{N_SEQS}", file=sys.stderr)
    print(f"Wrote {N_SEQS} sequences to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
