"""Exp 027 — cCRE balanced over SPECIALIZED classes only (5 classes, 10k each).

Exp 024 (8-class balanced) = 0.6921; Exp 026 (7-class no dELS) = 0.6925 (flat).
That rules out dELS as contributing distinct grammar (when its share is ≤12.5%).

This experiment drops the 3 "broad" enhancer/accessible classes (dELS, pELS, CA)
keeping only 5 SPECIALIZED classes: PLS (promoters), CA-CTCF (CTCF anchors),
TF (TF anchors), CA-H3K4me3 (H3K4me3-marked accessible), CA-TF (TF+accessible
no H3K27ac). Each gets 10k.

Hypothesis (aggressive): if all 3 enhancer-like classes are mutually redundant
AND the lift in exp 024 was specifically from the specialized classes, then
specialized-only might LIFT.

Counter-hypothesis: enhancer grammar (the dominant regulatory mode in the
genome) is essential context; removing it entirely will HURT and tell us
specialization is complementary to enhancers, not a replacement.

Either result is informative.
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
SEED = 27
KEEP = {"PLS", "CA-CTCF", "TF", "CA-H3K4me3", "CA-TF"}

rng = np.random.default_rng(SEED)

print("loading cCREs...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
df = df[df["class"].isin(KEEP)].reset_index(drop=True)
print(f"  kept {KEEP}: {len(df):,}")
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
print(f"  per_class quota: {per_class} (+{remainder} extra)")

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
