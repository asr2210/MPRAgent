"""Exp 008 — DHS topic-balanced across 12 non-label-aligned topics.

Forces equal representation per topic (~4167 per topic for 12 topics).
Tests whether explicit topic-balancing on top of exclusion further improves
diversity-driven gains.

Comparisons:
  - exp 005 (uniform from 12 topics, natural topic frequencies): 0.6752
  - baseline `dhs_stratified` (equal per topic, all 16 topics): 0.7055
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
SEED = 8

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)

topics = df["component"].unique()
n_per_topic = N // len(topics)
print(f"  topics: {len(topics)}, target per topic: {n_per_topic}")

pool = []
for topic in topics:
    sub = df[df["component"] == topic]
    print(f"  {topic}: {len(sub):,} avail")
    need = n_per_topic
    got = []
    attempt = 0
    while len(got) < need and attempt < 6:
        rem = need - len(got)
        idx = rng.choice(len(sub), size=rem * 2, replace=False if rem * 2 <= len(sub) else True)
        for i in idx:
            if len(got) >= need:
                break
            r = sub.iloc[i]
            seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
            if "N" in seq or len(seq) != L:
                continue
            got.append(seq)
        attempt += 1
    pool.extend(got)

# Top up to N with uniform from remaining if short
while len(pool) < N:
    r = df.iloc[rng.integers(0, len(df))]
    seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
    if "N" not in seq and len(seq) == L:
        pool.append(seq)
pool = pool[:N]

# Shuffle so training sees balanced mix
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {len(pool)} sequences x {L}bp to {OUT}")
