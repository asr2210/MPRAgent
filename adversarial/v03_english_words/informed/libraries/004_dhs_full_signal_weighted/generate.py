"""
Experiment 004: Full DHS Index, sampled with weight ∝ mean_signal.

Switch from the narrow 160k synthseqs subset to the full 3.6M-element DHS
Index. Sample with probability proportional to `mean_signal` (DHS signal
intensity averaged across biosamples) — strong-signal elements are more
likely to be active in our three measured cell types and contribute real
learning signal.

Sequences are extracted summit-centered, 200bp, from hg38. Elements where
the 200bp window would include any N base are filtered out (very rare).
Tissue distribution: natural component proportions from the index.

Comparison: random_uniform = 0.42 in MY pipeline. Goal: beat that with a real
DHS pool. Successful if eval_01 ≥ 0.48 (≈15% relative improvement).
"""

import sys
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


def extract_centered(fa: Fasta, chrom: str, summit: int) -> str | None:
    start = summit - LEN // 2
    end = start + LEN
    if start < 0:
        return None
    try:
        seq = fa[chrom][start:end]
    except (KeyError, IndexError):
        return None
    if len(seq) != LEN:
        return None
    seq = seq.upper()
    if "N" in seq:
        return None
    if not set(seq).issubset(set("ACGT")):
        return None
    return seq


def main() -> None:
    rng = np.random.default_rng(SEED)

    print("Loading DHS index...")
    df = pd.read_csv(DHS_TSV, sep="\t", low_memory=False)
    print(f"  {len(df):,} elements")

    # Restrict to standard chromosomes
    main_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    df = df[df["seqname"].isin(main_chroms)].reset_index(drop=True)
    print(f"  {len(df):,} after restricting to chr1-22,X,Y")

    # Build sampling weights ∝ mean_signal
    w = df["mean_signal"].astype(float).values
    p = w / w.sum()

    print("Loading hg38 fasta (pyfaidx, lazy)...")
    fa = Fasta(str(FASTA), as_raw=True, sequence_always_upper=True)

    # Oversample by 30% to cover any N-filtered or out-of-bounds losses
    target = N_TOTAL
    over = int(N_TOTAL * 1.3)
    print(f"Sampling {over} candidates (target {target})...")
    cand_idx = rng.choice(len(df), size=over, replace=False, p=p)
    cand = df.iloc[cand_idx]

    seqs = []
    components_used = []
    for _, row in cand.iterrows():
        s = extract_centered(fa, row["seqname"], int(row["summit"]))
        if s is None:
            continue
        seqs.append(s)
        components_used.append(row["component"])
        if len(seqs) >= target:
            break

    if len(seqs) < target:
        # Top up with random uniform; should rarely trigger
        print(f"WARNING: only got {len(seqs)} natural sequences; topping up with random")
        bases = np.array(list("ACGT"))
        while len(seqs) < target:
            arr = rng.integers(0, 4, size=LEN)
            seqs.append("".join(bases[arr]))

    assert len(seqs) == target
    for s in seqs:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")

    # Diagnostics
    used = pd.Series(components_used)
    print("\nComponent distribution in final library:")
    print(used.value_counts())


if __name__ == "__main__":
    main()
