"""Exp 007 — DHS from non-label-aligned topics, weighted by `numsamples`.

`numsamples` = number of biosamples (of 733) a DHS is detected in. Weighting by
this concentrates the library on broadly accessible (biologically robust)
regulatory elements — many tissue-invariant + multi-tissue shared.

Tests whether biologically robust DHSs train a better model than uniformly
sampled DHSs (which are dominated by single-biosample peaks).

Distribution of numsamples in full pool:
  1 sample:    ~1.23M DHSs
  2 samples:   ~430k
  3 samples:   ~250k
  ...
  733 samples: ~600 (tissue-invariant)
Weighted sampling shifts mass toward broadly accessible DHSs.

Predicts: eval_01 lifts above exp 005's 0.6752 if robustness helps; drops if
diversity-of-cell-type-specific is what matters.
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
SEED = 7

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index (with numsamples)...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component", "numsamples"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
print(f"  pool size: {len(df):,}")

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)

w = df["numsamples"].to_numpy().astype(float)
w = w / w.sum()
print(f"  weight mean={w.mean():.4e}, max={w.max():.4e}, min={w.min():.4e}")
print(f"  expected mean numsamples in sample: {(df['numsamples']*w).sum():.1f}")

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(df), size=need * 2, p=w, replace=True)
    for i in idx:
        if len(pool) >= N:
            break
        r = df.iloc[i]
        seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
        if "N" in seq or len(seq) != L:
            continue
        pool.append(seq)
    attempt += 1
assert len(pool) == N

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
