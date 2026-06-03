"""
Experiment 001: DHS sampling, signal-weighted, summit-centered.

Hypothesis: a strong informed baseline. DHS regions weighted by mean DNase
signal (with mild boost for biosample breadth) and centered on the summit
provides high-quality regulatory sequence coverage.

This is intended to roughly reproduce the dhs_topic baseline (~0.72 eval_01),
calibrating my pipeline before iterating.
"""
import gzip
import os
import numpy as np
import twobitreader

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DHS_FILE = os.path.join(ROOT, "data", "DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
TWOBIT = os.path.join(ROOT, "data", "hg38.2bit")
OUT = os.path.join(HERE, "sequences_0.txt")

WINDOW = 200
N_SEQ = 50_000
SEED = 1

# Standard primary assembly chromosomes only
CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_dhs():
    rows = []
    with gzip.open(DHS_FILE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[col["seqname"]] not in CHROMS:
                continue
            rows.append(
                (
                    p[col["seqname"]],
                    int(p[col["summit"]]),
                    float(p[col["mean_signal"]]),
                    int(p[col["numsamples"]]),
                    p[col["component"]],
                )
            )
    return rows


def sample_indices(weights, n, seed):
    rng = np.random.default_rng(seed)
    p = weights / weights.sum()
    return rng.choice(len(p), size=n, replace=False, p=p)


def extract(rows, idxs, genome):
    half = WINDOW // 2
    out = []
    for i in idxs:
        chrom, summit, _, _, _ = rows[i]
        start = max(0, summit - half)
        end = start + WINDOW
        seq = genome[chrom][start:end].upper()
        if len(seq) != WINDOW or any(c not in "ACGT" for c in seq):
            out.append(None)
        else:
            out.append(seq)
    return out


def main():
    print(f"Loading DHS index from {DHS_FILE}")
    rows = load_dhs()
    print(f"Loaded {len(rows):,} DHS rows on primary chromosomes")

    # Weight: mean_signal * log(numsamples + 1).
    # Strong signal AND mild boost for breadth across biosamples.
    weights = np.array(
        [r[2] * np.log(r[3] + 1) for r in rows], dtype=np.float64
    )

    print(f"Sampling {N_SEQ} with weighted DHS (seed={SEED})")
    # Oversample by 20% to allow rejecting N-containing windows.
    idxs = sample_indices(weights, int(N_SEQ * 1.2), seed=SEED)

    print("Opening 2bit genome and extracting sequences")
    genome = twobitreader.TwoBitFile(TWOBIT)
    seqs = extract(rows, idxs, genome)

    accepted = [s for s in seqs if s is not None][:N_SEQ]
    print(f"Accepted {len(accepted)} sequences")

    # If not enough (extremely unlikely), draw more
    if len(accepted) < N_SEQ:
        extra_idx = sample_indices(weights, N_SEQ, seed=SEED + 999)
        extra = extract(rows, extra_idx, genome)
        accepted += [s for s in extra if s is not None]
        accepted = accepted[:N_SEQ]

    assert len(accepted) == N_SEQ, f"Got {len(accepted)} not {N_SEQ}"
    assert all(len(s) == WINDOW for s in accepted)

    with open(OUT, "w") as f:
        for s in accepted:
            f.write(s + "\n")
    print(f"Wrote {len(accepted)} sequences to {OUT}")


if __name__ == "__main__":
    main()
