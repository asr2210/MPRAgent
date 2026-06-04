"""Exp 020 — DHS sampled proportional to MAX NMF loading per DHS.

Alternative interpretation of `dhs_topic` baseline: instead of summing all
16 topic loadings per DHS, weight by the maximum single-topic loading.
This more directly upweights "elements with strong cell-type-specific
accessibility signal" — a DHS with one very high loading (specific to one
topic) gets more weight than a DHS with diffuse loadings.

If exp 020 > exp 019 (sum-weighted, 0.6481) → max interpretation is closer
to baseline. If both < 0.6752 (uniform exp 005), the gap is elsewhere.
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
SEED = 20

rng = np.random.default_rng(SEED)

print("loading DHS index + NMF loadings...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
m = np.load(DATA / "NMF_Mixture.npy")
w = m.max(axis=0)  # per-DHS max loading
print(f"  weight (max loading) stats: mean={w.mean():.4f}, max={w.max():.4f}")

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
ok_arr = ok.to_numpy()

w_ok = np.where(ok_arr, w, 0.0)
p = w_ok / w_ok.sum()

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
