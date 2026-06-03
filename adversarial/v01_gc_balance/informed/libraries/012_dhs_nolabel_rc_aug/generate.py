"""Exp 012 — reverse-complement augmentation per DHS source.

H5 implies: don't concentrate / re-weight; augment per-source instead.

Method: sample 25k DHSs uniformly from the non-label-aligned pool (same recipe
as exp 005). For each DHS, emit BOTH the forward window and its reverse
complement. Library = 50k sequences from 25k biological sources × 2 strands.

If the model already learns RC-equivariance internally (likely for conv nets
with RC layers) → flat result.
If it doesn't → free lift from doubled effective strand coverage.

Comparison: exp 005 (50k DHSs, single strand) eval_01 = 0.6752.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"

N, L = 50_000, 200
HALF = L // 2
N_DHS = N // 2  # 25k DHSs, each gives 2 sequences (fwd + RC)
SEED = 12

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  pool size: {len(df):,}")

COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(COMP)[::-1]

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need_dhs = (N - len(pool) + 1) // 2
    idx = rng.choice(len(df), size=need_dhs * 2, replace=False if need_dhs * 2 <= len(df) else True)
    for i in idx:
        if len(pool) >= N:
            break
        r = df.iloc[i]
        seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
        if "N" in seq or len(seq) != L:
            continue
        pool.append(seq)
        if len(pool) < N:
            pool.append(revcomp(seq))
    attempt += 1
assert len(pool) == N

rng.shuffle(pool)  # interleave fwd/RC so training sees mixed batches

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
