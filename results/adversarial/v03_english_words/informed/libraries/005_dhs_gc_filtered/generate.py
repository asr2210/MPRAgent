"""
Experiment 005: DHS sequences filtered to GC 45-55%.

Tests the hypothesis that the underperformance of natural DHS sequences vs
random uniform in MY pipeline is dominated by GC variance. Random uniform has
tight ~50% GC; full-DHS has 47% ± 10%; baseline gc_50 (50% GC) beats both
random_uniform and gc_sweep — suggesting tight 50% GC is favored by this
pipeline.

If GC is the issue: a GC-filtered DHS library should beat both random_uniform
(0.42) AND full DHS (0.39), reaching maybe 0.43-0.46.

If GC is NOT the issue: result will stay near 0.39, suggesting deeper bias
(sequence composition beyond GC).
"""

from pathlib import Path
import numpy as np
import pandas as pd
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = Path(__file__).parent / "sequences_0.txt"
DHS_TSV = DATA / "DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz"
FASTA = DATA / "hg38.fa"

N_TOTAL = 50_000
LEN = 200
SEED = 0
GC_LOW = 0.45
GC_HIGH = 0.55


def extract_centered(fa, chrom, summit):
    start = summit - LEN // 2
    end = start + LEN
    if start < 0:
        return None
    try:
        seq = fa[chrom][start:end]
    except (KeyError, IndexError):
        return None
    if len(seq) != LEN or "N" in seq or not set(seq).issubset(set("ACGT")):
        return None
    return seq


def gc(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s)


def main() -> None:
    rng = np.random.default_rng(SEED)
    print("Loading DHS index...")
    df = pd.read_csv(DHS_TSV, sep="\t", low_memory=False)
    main_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    df = df[df["seqname"].isin(main_chroms)].reset_index(drop=True)
    print(f"  {len(df):,} after chrom filter")

    # Sample with mean_signal weighting; oversample heavily to absorb GC filter loss
    # We expect ~25-30% to pass the 45-55% GC filter, so oversample 4x
    over = N_TOTAL * 5
    w = df["mean_signal"].astype(float).values
    p = w / w.sum()
    print(f"Sampling {over} candidates...")
    cand_idx = rng.choice(len(df), size=min(over, len(df)), replace=False, p=p)
    cand = df.iloc[cand_idx]

    print("Loading hg38...")
    fa = Fasta(str(FASTA), as_raw=True, sequence_always_upper=True)

    seqs = []
    n_scanned = 0
    for _, row in cand.iterrows():
        n_scanned += 1
        s = extract_centered(fa, row["seqname"], int(row["summit"]))
        if s is None:
            continue
        if not (GC_LOW <= gc(s) <= GC_HIGH):
            continue
        seqs.append(s)
        if len(seqs) >= N_TOTAL:
            break

    print(f"After scanning {n_scanned}, got {len(seqs)} GC-filtered sequences")
    if len(seqs) < N_TOTAL:
        raise SystemExit(f"Not enough sequences: got {len(seqs)}, need {N_TOTAL}. "
                         "Increase oversample or relax GC bounds.")

    for s in seqs:
        assert len(s) == LEN and set(s).issubset(set("ACGT"))
    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")

    # Diagnostics
    gcs = np.array([gc(s) for s in seqs])
    print(f"GC: mean={gcs.mean():.3f} std={gcs.std():.3f} "
          f"min={gcs.min():.3f} max={gcs.max():.3f}")


if __name__ == "__main__":
    main()
