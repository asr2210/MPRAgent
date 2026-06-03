"""
007_composite_mix
=================
Maximum-diversity composite library:
  - 12,500 sequences: cCRE-derived (200bp around midpoint)
  - 12,500 sequences: DHS specificity-weighted (200bp around summit)
  - 12,500 sequences: Random hg38 genomic windows (no annotation)
  - 12,500 sequences: Synthetic random + 4-8 planted JASPAR motifs

Different evals appear to favor different library types. This composite
maximizes coverage across all four "axes" probed so far.

Hypothesis: if eval is a heterogeneous mix that benefits from broad source
coverage, this scores ~max of components (~0.05). If different libraries
contribute orthogonal signal, this might score higher than any single
component.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
QUARTER = 12_500
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CCRE_PATH = DATA_DIR / "ccre_v4.bed"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
JASPAR_PATH = DATA_DIR / "jaspar2024_core.txt"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_ccre_recs():
    recs = []
    with open(CCRE_PATH) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            mid = (int(parts[1]) + int(parts[2])) // 2
            recs.append((parts[0], mid))
    return recs


def load_dhs_recs():
    recs = []
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            recs.append((parts[0], int(parts[6]), int(parts[5])))
    return recs


def parse_jaspar(path):
    pfms = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith(">"):
            rows = []
            for j in range(1, 5):
                row = lines[i + j].strip()
                start = row.index("[") + 1
                end = row.index("]")
                rows.append([float(x) for x in row[start:end].split()])
            mat = np.array(rows, dtype=np.float64)
            col_sums = mat.sum(axis=0, keepdims=True)
            col_sums[col_sums == 0] = 1.0
            mat = (mat + 0.01) / (col_sums + 0.04)
            if 6 <= mat.shape[1] <= 20:
                pfms.append(mat)
            i += 5
        else:
            i += 1
    return pfms


def extract_window(fa, chrom, center, half, seq_len):
    L = len(fa[chrom])
    start = max(0, center - half)
    end = start + seq_len
    if end > L:
        end = L
        start = end - seq_len
    if start < 0:
        return None
    s = str(fa[chrom][start:end])
    if len(s) != seq_len or set(s) - set("ACGT"):
        return None
    return s


def sample_window_seqs(recs, n, fa, rng, key="center"):
    seqs = []
    perm = rng.permutation(len(recs))
    pi = 0
    while len(seqs) < n and pi < len(perm):
        rec = recs[perm[pi]]
        pi += 1
        chrom = rec[0]
        center = rec[1]
        s = extract_window(fa, chrom, center, HALF, SEQ_LEN)
        if s:
            seqs.append(s)
    return seqs


def sample_random_genomic(n, fa, rng):
    chr_names = list(CANONICAL)
    chr_lens = {c: len(fa[c]) for c in chr_names}
    probs = np.array([chr_lens[c] for c in chr_names], dtype=np.float64)
    probs /= probs.sum()
    seqs = []
    while len(seqs) < n:
        batch = max(500, n - len(seqs))
        chr_idx = rng.choice(len(chr_names), size=batch, p=probs)
        for i in chr_idx:
            chrom = chr_names[i]
            L = chr_lens[chrom]
            start = int(rng.integers(0, L - SEQ_LEN + 1))
            s = str(fa[chrom][start:start + SEQ_LEN])
            if len(s) != SEQ_LEN or set(s) - set("ACGT"):
                continue
            seqs.append(s)
            if len(seqs) >= n:
                break
    return seqs


def sample_motif_synth(n, pfms, rng):
    bases = np.array(["A", "C", "G", "T"])
    seqs = []
    for _ in range(n):
        bg = list(rng.choice(bases, size=SEQ_LEN))
        n_motifs = int(rng.integers(4, 9))
        chosen = rng.choice(len(pfms), size=n_motifs, replace=True)
        occ = []
        for k in chosen:
            pfm = pfms[k]
            ml = pfm.shape[1]
            inst = "".join(rng.choice(bases, p=pfm[:, c]) for c in range(ml))
            for _try in range(40):
                st = int(rng.integers(0, SEQ_LEN - ml + 1))
                en = st + ml
                if all(en <= a or st >= b for a, b in occ):
                    break
            else:
                st = int(rng.integers(0, SEQ_LEN - ml + 1))
                en = st + ml
            occ.append((st, en))
            for off, ch in enumerate(inst):
                bg[st + off] = ch
        seqs.append("".join(bg))
    return seqs


def main():
    rng = np.random.default_rng(SEED)
    print("Loading data sources...", file=sys.stderr)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)
    ccres = load_ccre_recs()
    dhs = load_dhs_recs()
    # Reweight DHS by specificity
    nsamples = np.array([r[2] for r in dhs])
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    # Pre-pick DHS indices via weighted sampling
    pfms = parse_jaspar(JASPAR_PATH)
    print(f"  ccres={len(ccres)}, dhs={len(dhs)}, pfms={len(pfms)}", file=sys.stderr)

    print("Sampling 12500 cCRE windows...", file=sys.stderr)
    ccre_seqs = sample_window_seqs(ccres, QUARTER, fa, rng)
    print("Sampling 12500 DHS specificity-weighted windows...", file=sys.stderr)
    dhs_pool = rng.choice(len(dhs), size=int(QUARTER * 1.3), replace=False, p=w)
    dhs_seqs = []
    for i in dhs_pool:
        if len(dhs_seqs) >= QUARTER:
            break
        chrom, summit, _ = dhs[i]
        s = extract_window(fa, chrom, summit, HALF, SEQ_LEN)
        if s:
            dhs_seqs.append(s)
    print(f"  got {len(dhs_seqs)} DHS", file=sys.stderr)

    print("Sampling 12500 random genomic windows...", file=sys.stderr)
    rg_seqs = sample_random_genomic(QUARTER, fa, rng)
    print("Generating 12500 motif synthetic sequences...", file=sys.stderr)
    motif_seqs = sample_motif_synth(QUARTER, pfms, rng)

    all_seqs = ccre_seqs + dhs_seqs + rg_seqs + motif_seqs
    assert len(all_seqs) == N_SEQS, f"got {len(all_seqs)}"
    rng.shuffle(all_seqs)

    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in all_seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
