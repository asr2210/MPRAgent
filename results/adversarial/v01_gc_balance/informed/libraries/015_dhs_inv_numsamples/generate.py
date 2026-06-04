"""Exp 015 — DHS weighted INVERSELY by `numsamples` (favor cell-type-specific).

Re-read baseline recipes: `dhs_topic` (0.7232, best baseline) samples
proportional to NMF topic loadings, which UPWEIGHTS cell-type-specific DHSs
(those with strong topic signal). Exp 007 weighted by `numsamples` (the
OPPOSITE direction — favoring broadly-accessible DHSs) and failed (0.6185).

H5 ("concentration always hurts") was wrong: direction matters. Now H6:
> Concentration on the cell-type-SPECIFIC axis helps. Concentration on the
> cell-type-INVARIANT axis hurts.

Test: weight each DHS by 1/(1 + numsamples). DHSs detected in few biosamples
(low numsamples = cell-type-specific) get high weight; broadly-accessible
DHSs get tiny weight. Combine with the 12-topic exclusion (best filter so far).

Predicts: eval_01 > 0.6752. If it lifts toward the 0.7232 baseline, H6 is
strongly supported — direction matters and inverse weighting is the right move.
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
SEED = 15

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index (with numsamples)...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component", "numsamples"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
print(f"  pool size: {len(df):,}")
print(f"  numsamples distribution: median={int(df['numsamples'].median())}, "
      f"p90={int(df['numsamples'].quantile(0.9))}, max={int(df['numsamples'].max())}")

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)

w = 1.0 / (1.0 + df["numsamples"].to_numpy().astype(float))
w = w / w.sum()
exp_num = (df["numsamples"] * w).sum()
print(f"  expected mean numsamples after inverse-weight: {exp_num:.2f} "
      f"(uniform would be ~{df['numsamples'].mean():.2f})")

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(df), size=need * 2, p=w, replace=True)
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
