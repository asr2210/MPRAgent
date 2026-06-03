"""Exp 029 — cCRE 7-class (no dELS) × 23-chrom balanced.

Combines wins from:
- exp 026 (drop dELS, 7-class balanced) = 0.6925 (tied best, flat vs exp 024)
- exp 028 (8-class × chrom-balanced) = 0.6940 NEW BEST

Tests: does dropping dELS still produce a flat result when chrom-balance is
also applied? Or does chrom-balance need dELS's volume to maintain coverage?

7 classes × 23 chroms × ~310/cell = 50,030 (trimmed to 50k).
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
SEED = 29
DROP = {"dELS"}
USE_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

rng = np.random.default_rng(SEED)

print("loading cCREs...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
df = df[~df["class"].isin(DROP)].reset_index(drop=True)
df = df[df["chrom"].isin(USE_CHROMS)].reset_index(drop=True)
print(f"  after drop+chrom: {len(df):,}")

df["mid"] = (df["start"] + df["end"]) // 2
df["win_start"] = df["mid"] - HALF
df["win_end"] = df["mid"] + HALF

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["chrom"].unique()}
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
df = df[ok].reset_index(drop=True)

classes = sorted(df["class"].unique())
n_cells = len(classes) * len(USE_CHROMS)
per_cell = N // n_cells
remainder = N - per_cell * n_cells
print(f"  classes={classes}, per_cell={per_cell}, remainder={remainder}")

pool = []
extras = remainder
for cls in classes:
    sub_cls = df[df["class"] == cls]
    for chrom in USE_CHROMS:
        sub = sub_cls[sub_cls["chrom"] == chrom].reset_index(drop=True)
        quota = per_cell + (1 if extras > 0 else 0)
        if extras > 0:
            extras -= 1
        if len(sub) == 0:
            sub = sub_cls.reset_index(drop=True)
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
        if len(collected) < quota:
            sub_all = sub_cls.reset_index(drop=True)
            while len(collected) < quota:
                r = sub_all.iloc[rng.integers(len(sub_all))]
                seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
                if "N" in seq or len(seq) != L:
                    continue
                collected.append(seq)
        pool.extend(collected)

assert len(pool) == N, f"got {len(pool)}"
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
