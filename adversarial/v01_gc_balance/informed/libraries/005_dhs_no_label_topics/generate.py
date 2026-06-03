"""Exp 005 — DHS uniform, EXCLUDING topics aligned with labeling cell types.

Direct test of generalisation. The labeling cells (K562, HepG2, SK-N-SH) map to:
  - Myeloid/erythroid (K562)
  - Digestive (HepG2 — hepatic carcinoma)
  - Cancer/epithelial (HepG2 — secondary)
  - Neural (SK-N-SH)

Sample 50k uniformly from DHSs whose dominant topic is NONE of these four.
That leaves ~2.6M DHSs across 12 'other' regulatory programs.

Prediction (H2): if the model learns universal regulatory grammar, eval_01
should match or barely lag pure DHS uniform (~0.66). If labeling-cell-aligned
regulatory programs are needed, eval_01 will drop substantially.
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
SEED = 5

EXCLUDED_TOPICS = {
    "Myeloid / erythroid",   # K562
    "Digestive",             # HepG2 (hepatic)
    "Cancer / epithelial",   # HepG2 (carcinoma secondary)
    "Neural",                # SK-N-SH
}

print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
print(f"  total DHSs: {len(df):,}")
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
print(f"  after excluding label-aligned topics: {len(df):,}")
print("  remaining topic counts:")
print(df["component"].value_counts())

print("loading hg38...")
fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}

df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)

rng = np.random.default_rng(SEED)
pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(df), size=need * 2, replace=False if need * 2 <= len(df) else True)
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
