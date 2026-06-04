"""Exp 026 — cCRE balanced EXCLUDING dELS (7 classes, ~7143 each).

Exp 024 (cCRE 8-class balanced, dELS at 12.5%) = 0.6921 NEW BEST.
Exp 022 (cCRE all-classes uniform, dELS at 62.5%) = 0.6827.

The lift from 022 → 024 came from REDISTRIBUTING toward rare classes
(PLS, CA-TF, CA-H3K4me3, TF, CA-CTCF). dELS dropped from 62.5% to 12.5%.

Question: does dELS still contribute, or is its grammar already covered by
pELS (which is also "enhancer-like")? Test by dropping dELS entirely and
balancing the remaining 7 classes: CA, CA-CTCF, CA-H3K4me3, CA-TF, PLS,
TF, pELS at ~7143 each.

- If eval_01 > 0.6921 → dELS was still net redundant; rare-class is even
  better at higher share.
- If eval_01 ≈ 0.6921 → dELS at 12.5% adds nothing on top of the others.
- If eval_01 < 0.6921 → dELS does contribute distinct grammar; drop hurts.
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
SEED = 26
DROP = {"dELS"}

rng = np.random.default_rng(SEED)

print("loading cCREs...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
df = df[~df["class"].isin(DROP)].reset_index(drop=True)
print(f"  after dropping {DROP}: {len(df):,}")
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

classes = sorted(df["class"].unique())
per_class = N // len(classes)
remainder = N - per_class * len(classes)
print(f"  classes: {classes}")
print(f"  per_class quota: {per_class} (+{remainder} extra to first)")

pool = []
for ci, cls in enumerate(classes):
    quota = per_class + (1 if ci < remainder else 0)
    sub = df[df["class"] == cls].reset_index(drop=True)
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
    assert len(collected) == quota
    pool.extend(collected)
    print(f"  {cls}: pool {len(sub):,}, took {quota}")

assert len(pool) == N
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
