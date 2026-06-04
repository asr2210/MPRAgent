"""Exp 024 — cCRE class-balanced (equal samples per of 8 classes).

Exp 022 (cCRE all-classes uniform) = 0.6827 sampled from the natural class
distribution: dELS 62.5%, pELS 10.6%, CA 10.5%, CA-CTCF 5.4%, TF 4.5%,
CA-H3K4me3 3.4%, PLS 2.0%, CA-TF 1.1%. With 50k draws this means dELS
dominates (~31k of 50k) and PLS+CA-TF together only ~1.5k.

This experiment forces equal representation: ~6250 sequences per class.
Tests whether natural proportions are optimal or whether the rare classes
(PLS promoters, CA-TF) carry distinctive regulatory grammar that we're
under-sampling.

Predicts: if rare classes hold unique signal → lift over exp 022. If dELS
density was already the right mix → flat or hurt.
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
SEED = 24

rng = np.random.default_rng(SEED)

print("loading cCREs (all classes)...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
print(f"  total: {len(df):,}")
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
print(f"  after edge filter: {len(df):,}")

classes = sorted(df["class"].unique())
per_class = N // len(classes)
remainder = N - per_class * len(classes)
print(f"  classes: {classes}")
print(f"  per_class quota: {per_class} (+{remainder} distributed)")

pool = []
for ci, cls in enumerate(classes):
    quota = per_class + (1 if ci < remainder else 0)
    sub = df[df["class"] == cls].reset_index(drop=True)
    print(f"  {cls}: pool {len(sub):,}, draw {quota}")
    collected = []
    attempt = 0
    while len(collected) < quota and attempt < 10:
        need = quota - len(collected)
        size = min(len(sub), need * 3)
        replace = need * 3 > len(sub)
        idx = rng.choice(len(sub), size=size, replace=replace)
        for i in idx:
            if len(collected) >= quota:
                break
            r = sub.iloc[i]
            seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
            if "N" in seq or len(seq) != L:
                continue
            collected.append(seq)
        attempt += 1
    assert len(collected) == quota, f"{cls}: got {len(collected)}/{quota}"
    pool.extend(collected)

assert len(pool) == N
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
