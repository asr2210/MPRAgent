"""Exp 014 — looser topic exclusion: exclude only the 2 most label-specific topics.

Topic-to-label mapping:
  Myeloid / erythroid  ~  K562 (chronic myelogenous leukemia / erythroleukemia)
  Digestive            ~  HepG2 (hepatocellular carcinoma, liver = digestive)
  Cancer / epithelial  ~  HepG2 (HCC also has cancer-epithelial signature)
  Neural               ~  SK-N-SH (neuroblastoma)

Exp 005 removed all 4 (kept 12 topics). This experiment removes only the 2
MOST specific (Myeloid/erythroid + Neural). 14 topics retained.

Question: is the lift from exp 005 (vs exp 003 all 16) driven by removing
ALL label-aligned, or is removing just the two clearly-specific enough?

If exp 014 ≥ exp 005 (0.6752) → looser exclusion sufficient; redundant cancer/digestive removal was unnecessary
If exp 014 < exp 005 → 4-topic exclusion does the work — all label-aligned regulatory programs need to leave
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
SEED = 14

EXCLUDED_TOPICS = {"Myeloid / erythroid", "Neural"}

rng = np.random.default_rng(SEED)

print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
print(f"  pool: {len(df):,} (14 topics)")

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)

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
