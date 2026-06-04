"""Exp 021 — combine no-label-topic filter + NMF max-loading weighting.

Best filter (exp 005, 0.6752): exclude 4 label-aligned topics, uniform within.
Best NMF interpretation (exp 020, 0.6641): max loading; underperforms uniform.

This combines: same filter (12 non-label topics) + within that pool, weight
by per-DHS MAX NMF loading (cell-specificity score). Tests whether the
combination unlocks something neither alone does.

If eval_01 > 0.6752 → combination helps; cell-specificity weighting works
when the labeling-aligned topics are first removed.
If <= 0.6752 → uniform-within-filter remains the local optimum.
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
SEED = 21

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index + NMF loadings...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
m = np.load(DATA / "NMF_Mixture.npy")
w_all = m.max(axis=0)
print(f"  total DHSs: {len(df):,}; matrix cols: {m.shape[1]}")

# Filter by topic
keep_mask = (~df["component"].isin(EXCLUDED_TOPICS)).to_numpy()
print(f"  no-label-aligned mask: {keep_mask.sum():,}")

# Edge filter
fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
ok_arr = ok.to_numpy()

mask = keep_mask & ok_arr
print(f"  pool after topic + edge filter: {mask.sum():,}")

w = np.where(mask, w_all, 0.0)
p = w / w.sum()

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(df), size=need * 2, p=p, replace=True)
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
