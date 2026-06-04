"""Exp 028 — cCRE 8-class balanced + per-chromosome balanced within class.

Exp 024 (8-class balanced, uniform within class) = 0.6921.
Within each class, sampling is uniform over the full pool — but the genome
is highly non-uniform in cCRE density per chromosome. chr1 has ~10x more
cCREs than chr21. So uniform draws within class effectively oversample
chr1/chr2/chr3 etc and undersample small chroms.

This experiment adds a second orthogonal balance: within each class, draw
the SAME number from each chromosome (or as close as possible). 22 autosomes
+ X = 23 chroms × 8 classes × ~272/cell ≈ 50k.

If chrom variation in cCRE density correlates with biology (e.g. tissue-specific
hotspots) → chrom-balance hurts because we're flattening real signal.
If chrom variation is mostly genomic-length confounding → chrom-balance helps.
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
SEED = 28

rng = np.random.default_rng(SEED)

print("loading cCREs...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
df["mid"] = (df["start"] + df["end"]) // 2
df["win_start"] = df["mid"] - HALF
df["win_end"] = df["mid"] + HALF

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chroms_fa = set(fa.keys())
df = df[df["chrom"].isin(chroms_fa)].reset_index(drop=True)
chrom_len = {c: len(fa[c]) for c in df["chrom"].unique()}
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
df = df[ok].reset_index(drop=True)

# Only autosomes + X (drop Y/M; very sparse cCRE coverage)
USE_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
df = df[df["chrom"].isin(USE_CHROMS)].reset_index(drop=True)

classes = sorted(df["class"].unique())
n_classes = len(classes)
n_chroms = len(USE_CHROMS)
per_cell = N // (n_classes * n_chroms)
remainder = N - per_cell * n_classes * n_chroms
print(f"  per_cell base quota: {per_cell} (+{remainder} extras)")

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
            # No cCREs of this class on this chrom — borrow from any chrom for class
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
            # Fill from class pool
            sub_all = sub_cls.reset_index(drop=True)
            while len(collected) < quota:
                r = sub_all.iloc[rng.integers(len(sub_all))]
                seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
                if "N" in seq or len(seq) != L:
                    continue
                collected.append(seq)
        pool.extend(collected)

assert len(pool) == N, f"got {len(pool)} expected {N}"
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
