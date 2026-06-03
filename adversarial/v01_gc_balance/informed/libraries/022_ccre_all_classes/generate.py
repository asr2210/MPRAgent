"""Exp 022 — ENCODE cCRE all-classes uniform.

Exp 010 used dELS-only cCREs (1.47M) and tied with exp 005. This experiment
uses ALL cCRE classes (dELS, pELS, PLS, CA-CTCF, TF, CA-H3K4me3, CA, CA-TF
= 2.35M total). Tests whether including promoters, CTCF sites, TF-only, and
chromatin-accessible elements alongside distal enhancers adds class diversity
that helps.

Predicts: comparable to or slightly better than exp 010 if class diversity
helps; flat or worse if pELS/PLS/CTCF add nothing for MPRA prediction.
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
SEED = 22

rng = np.random.default_rng(SEED)

print("loading cCREs (all classes)...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
print(f"  total: {len(df):,}; class distribution:")
print(df["class"].value_counts().to_string())

df["mid"] = (df["start"] + df["end"]) // 2
df["win_start"] = df["mid"] - HALF
df["win_end"] = df["mid"] + HALF

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chroms = set(fa.keys())
df = df[df["chrom"].isin(chroms)].reset_index(drop=True)
chrom_len = {c: len(fa[c]) for c in df["chrom"].unique()}
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
df = df[ok].reset_index(drop=True)

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
