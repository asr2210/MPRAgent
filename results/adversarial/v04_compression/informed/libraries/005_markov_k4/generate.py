#!/usr/bin/env python3
"""
005_markov_k4 — synthetic sequences with hg38 4-mer composition.

50k 200bp sequences generated from a 4th-order Markov chain trained on
hg38 chr1, chr17, chr19, chr22 (the same chromosomes used in 002 for
the real-DNA baseline).

Purpose: decompose what about real DNA gives it the +0.19 r-point
advantage over uniform synthetic. If matching 4-mer composition is
sufficient, this experiment will reproduce the real-DNA score. If
not, higher-order structure (longer motifs, repeats, long-range
correlations) is what matters.

Compare to:
- 001 synth random (uniform i.i.d.): eval_01 = 0.3068
- 002 real genome (chr1/17/19/22):   eval_01 = 0.4992
- 005 markov_k4 (this):              eval_01 = ?

Interpretation cases:
- 005 ≈ 0.31: k=4 not enough; need higher order.
- 005 ≈ 0.40: partial; some signal is in k=4, more in higher.
- 005 ≈ 0.50: 4-mer composition is the whole real-DNA prior. Means
  synthetic sequences are sufficient if k-mer matched.

Generalization argument: if k-mer composition explains the real-DNA
prior, then ANY cell type's evaluation should benefit from the same
prior (because k-mer composition is a property of DNA, not biology).
A k-mer-matched library is just as transferable as real DNA — and
could be more interesting because it lets us decouple the prior from
incidental functional bias toward our 3 cell types.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
K = 4  # context length; 4th-order Markov ⇒ 256 contexts × 4 outputs
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CHROMS = ["chr1", "chr17", "chr19", "chr22"]
ALPHABET = "ACGT"
BASE2INT = {b: i for i, b in enumerate(ALPHABET)}


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        first = f.readline()
        assert first.startswith(">")
        return "".join(line.strip() for line in f).upper()


def fit_markov(seqs, k: int):
    """Estimate P(base | preceding k bases) on contiguous runs of ACGT."""
    counts = np.ones((4**k, 4), dtype=np.float64)  # Laplace smoothing
    init_counts = np.ones(4**k, dtype=np.float64)  # for initial k-mer
    valid = set("ACGT")
    for seq in seqs:
        # split on any non-ACGT character → contiguous runs
        run = []
        for ch in seq:
            if ch in valid:
                run.append(ch)
            else:
                _accumulate(run, k, counts, init_counts)
                run = []
        _accumulate(run, k, counts, init_counts)
    # normalize
    probs = counts / counts.sum(axis=1, keepdims=True)
    init_probs = init_counts / init_counts.sum()
    return probs, init_probs


def _accumulate(run, k, counts, init_counts):
    if len(run) < k + 1:
        return
    idx = 0
    for ch in run[:k]:
        idx = idx * 4 + BASE2INT[ch]
    init_counts[idx] += 1
    for ch in run[k:]:
        nxt = BASE2INT[ch]
        counts[idx, nxt] += 1
        # shift context: drop highest, multiply by 4, add new
        idx = (idx % (4 ** (k - 1))) * 4 + nxt


def sample_markov(rng, probs, init_probs, k: int, L: int, n: int):
    """Generate n sequences of length L."""
    out = np.empty((n, L), dtype=np.int8)
    # vectorized initial k-mer sample
    init_ctx = rng.choice(len(init_probs), size=n, p=init_probs)
    # decode k initial bases
    for j in range(k):
        out[:, k - 1 - j] = (init_ctx // (4**j)) % 4
    ctx = init_ctx.copy()
    mask = 4 ** (k - 1)
    for j in range(k, L):
        # Sample next base for each row given ctx
        # Build cumulative for each row
        u = rng.random(n)
        row_probs = probs[ctx]  # (n, 4)
        cum = np.cumsum(row_probs, axis=1)
        nxt = (u[:, None] >= cum).sum(axis=1)
        out[:, j] = nxt
        ctx = (ctx % mask) * 4 + nxt
    return out


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for c in CHROMS:
        s = load_chrom(DATA_DIR / f"{c}.fa.gz")
        seqs.append(s)
        print(f"loaded {c}: {len(s):,} bp")
    print(f"fitting {K}-th order Markov chain...")
    probs, init_probs = fit_markov(seqs, K)
    print(f"  transitions shape {probs.shape}, init shape {init_probs.shape}")
    print(f"sampling {N} sequences of length {L}...")
    out = sample_markov(rng, probs, init_probs, K, L, N)
    alphabet_arr = np.array(list(ALPHABET))
    strs = ["".join(alphabet_arr[row]) for row in out]
    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(strs))
        f.write("\n")
    print(f"wrote {N} sequences → {out_path}")


if __name__ == "__main__":
    main()
