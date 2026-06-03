"""Exp 010 — ENCODE cCRE distal-enhancer-like-signature (dELS) library.

Orthogonal annotation source to DHS Index. ENCODE cCREs (Registry V4) are
curated cis-regulatory elements integrating DNase + H3K4me3 + H3K27ac + CTCF
across 1518 biosamples. Classes:
  dELS  1,469,205  distal enhancer-like signature
  pELS    249,464  proximal ELS
  CA      245,985  chromatin accessible only
  ...

This experiment: uniform sample 50k from dELS only — putative enhancers,
which is what MPRA primarily reports on. Tests if a curated enhancer-only
library beats the broader DHS pool (exp 005 = 0.6752).

cCRE widths are 150-350bp; we center on midpoint and take ±100bp window.
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
SEED = 10

rng = np.random.default_rng(SEED)

print("loading cCREs...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
df = df[df["class"] == "dELS"].reset_index(drop=True)
print(f"  dELS only: {len(df):,}")

df["mid"] = (df["start"] + df["end"]) // 2
df["win_start"] = df["mid"] - HALF
df["win_end"] = df["mid"] + HALF

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chroms = set(fa.keys())
df = df[df["chrom"].isin(chroms)].reset_index(drop=True)
chrom_len = {c: len(fa[c]) for c in df["chrom"].unique()}
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  after edge filter: {len(df):,}")

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(df), size=need * 2, replace=False if need * 2 <= len(df) else True)
    for i in idx:
        if len(pool) >= N:
            break
        r = df.iloc[i]
        seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
        if "N" in seq or len(seq) != L:
            continue
        pool.append(seq)
    attempt += 1
assert len(pool) == N

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
