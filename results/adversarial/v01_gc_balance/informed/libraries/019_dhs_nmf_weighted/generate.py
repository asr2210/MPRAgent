"""Exp 019 — DHS sampled proportional to NMF topic loadings (replicate baseline).

Downloaded the per-DHS NMF Mixture matrix from Meuleman Zenodo
(2018-06-08NC16_NNDSVD_Mixture.npy, shape 16 x 3,591,898, aligned 1-to-1 with
DHS Index rows).

`dhs_topic` baseline recipe: "sampled with probability proportional to NMF
topic loadings (16 topics), which upweights elements with strong cell-type-
specific accessibility signal." Simplest faithful reading: p_i ∝ sum of
loadings across all 16 topics per DHS.

This is the FIRST experiment with the full per-DHS topic information. If it
reproduces ~0.72 (baseline `dhs_topic` is 0.7232), we've closed the gap modulo
multi-seed averaging. If still ~0.67, the gap is fully attributable to seed
variance and there's nothing more to do at the recipe level.
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
SEED = 19

rng = np.random.default_rng(SEED)

print("loading DHS index + NMF loadings...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
m = np.load(DATA / "NMF_Mixture.npy")  # (16, N_dhs)
assert m.shape[1] == len(df), f"matrix cols {m.shape[1]} != index rows {len(df)}"
print(f"  DHSs: {len(df):,}; NMF matrix shape: {m.shape}")

w = m.sum(axis=0)  # per-DHS total loading
print(f"  weight stats: mean={w.mean():.4f}, std={w.std():.4f}, "
      f"p99={np.percentile(w, 99):.4f}, max={w.max():.4f}")

# Edge filter for window viability
fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
ok_arr = ok.to_numpy()
print(f"  edge-OK DHSs: {ok_arr.sum():,}")

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
