"""
013_dinuc_shuffle_dhs
=====================
50k DHS (specificity-weighted) sequences, each DINUCLEOTIDE-SHUFFLED.
Preserves: per-sequence GC content, dinucleotide frequencies.
Destroys: motifs, repeat structure, conservation patterns, all higher-order.

Probe: if eval ~= dhs_specific (0.049), eval is composition-driven and
indifferent to motif content. If eval crashes to << 0.04, eval depends
on real sequence motifs/structure.

Uses Altschul-Erickson dinucleotide shuffle via Eulerian random walk.
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def dinuc_shuffle(s: str, rng: np.random.Generator) -> str:
    """Altschul-Erickson dinucleotide-preserving shuffle via Eulerian walk."""
    n = len(s)
    if n < 2:
        return s
    # build edge list: for each base, the bases that follow it
    edges = defaultdict(list)
    for i in range(n - 1):
        edges[s[i]].append(s[i + 1])
    last = s[-1]
    # try up to 20 times to get a valid Eulerian path
    for _attempt in range(20):
        # shuffle edges per base (will yield random walks); reserve last
        # edge per source to land at `last`
        edge_shuf = {b: list(es) for b, es in edges.items()}
        # remove a final outgoing edge per source = last edge of last walk
        # Simpler: shuffle each list, then walk greedily; if dead-end mid-walk,
        # restart.
        for b in edge_shuf:
            rng.shuffle(edge_shuf[b])
        # Walk
        out = [s[0]]
        cur = s[0]
        ok = True
        for _ in range(n - 1):
            if not edge_shuf[cur]:
                ok = False
                break
            nxt = edge_shuf[cur].pop()
            out.append(nxt)
            cur = nxt
        if ok and len(out) == n:
            return "".join(out)
    return s  # fallback: return original (rare)


def load_dhs():
    chroms, summits, nsamples = [], [], []
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in CANONICAL:
                continue
            chroms.append(parts[0])
            summits.append(int(parts[6]))
            nsamples.append(int(parts[5]))
    return chroms, np.array(summits), np.array(nsamples, dtype=np.int32)


def main():
    rng = np.random.default_rng(SEED)
    chroms, summits, nsamples = load_dhs()
    w = 1.0 / np.sqrt(nsamples.astype(np.float64))
    w /= w.sum()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    pool = rng.choice(len(chroms), size=int(N_SEQS * 1.3), replace=False, p=w)
    seqs = []
    pi = 0
    while len(seqs) < N_SEQS and pi < len(pool):
        i = pool[pi]
        pi += 1
        chrom = chroms[i]
        summit = int(summits[i])
        L = len(fa[chrom])
        start = max(0, summit - HALF)
        end = start + SEQ_LEN
        if end > L:
            end = L
            start = end - SEQ_LEN
        if start < 0:
            continue
        s = str(fa[chrom][start:end])
        if len(s) != SEQ_LEN or set(s) - set("ACGT"):
            continue
        shuffled = dinuc_shuffle(s, rng)
        if len(shuffled) != SEQ_LEN:
            continue
        seqs.append(shuffled)
        if len(seqs) % 10_000 == 0:
            print(f"  {len(seqs)}/{N_SEQS}", file=sys.stderr)
    assert len(seqs) == N_SEQS
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
