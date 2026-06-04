"""Exp 023 — cCRE all-classes INTERSECTED with no-label-aligned DHS pool.

Exp 022 (cCRE all-classes uniform) = 0.6827 — new best.
Exp 005 (DHS no-label uniform) = 0.6752 — prior best.

This experiment combines both wins: take all cCREs, but only keep those that
overlap a DHS whose dominant NMF topic is NOT label-aligned (4-topic exclusion).
The intersection is high-confidence regulatory elements that don't carry
direct label-cell signal.

If eval_01 > 0.6827 → both filters help orthogonally.
If = 0.6827 → cCRE selection alone is enough; topic exclusion adds nothing
once you're on the cCRE catalog.
If < 0.6827 → topic exclusion HURTS cCREs (perhaps removes useful elements).

Implementation: per-chrom sort-merge to find which cCREs (by midpoint) fall
inside any no-label-aligned DHS interval.
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
SEED = 23

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index (no-label only)...")
dhs = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                  usecols=["seqname", "start", "end", "component"])
dhs = dhs[~dhs["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
print(f"  no-label DHSs: {len(dhs):,}")

print("loading cCREs (all classes)...")
ccre = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                   names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
ccre["mid"] = (ccre["start"] + ccre["end"]) // 2
print(f"  cCREs: {len(ccre):,}")

# Per-chrom interval test: is each cCRE midpoint inside any no-label DHS?
print("intersecting (per-chrom merge_asof)...")
keep_mask = np.zeros(len(ccre), dtype=bool)
for chrom in sorted(set(dhs["seqname"]).intersection(set(ccre["chrom"]))):
    d = dhs[dhs["seqname"] == chrom].sort_values("start").reset_index(drop=True)
    c_idx = np.where(ccre["chrom"] == chrom)[0]
    c_mid = ccre.iloc[c_idx]["mid"].to_numpy()
    # For each c_mid, find largest d.start <= c_mid
    d_starts = d["start"].to_numpy()
    d_ends = d["end"].to_numpy()
    pos = np.searchsorted(d_starts, c_mid, side="right") - 1
    valid = pos >= 0
    inside = np.zeros(len(c_mid), dtype=bool)
    inside[valid] = (d_ends[pos[valid]] >= c_mid[valid])
    keep_mask[c_idx] = inside
print(f"  cCREs whose midpoint is in a no-label DHS: {keep_mask.sum():,}")

ccre = ccre[keep_mask].reset_index(drop=True)
print(f"  class distribution after intersection:")
print(ccre["class"].value_counts().to_string())

ccre["win_start"] = ccre["mid"] - HALF
ccre["win_end"] = ccre["mid"] + HALF

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chroms = set(fa.keys())
ccre = ccre[ccre["chrom"].isin(chroms)].reset_index(drop=True)
chrom_len = {c: len(fa[c]) for c in ccre["chrom"].unique()}
ok = (ccre["win_start"] >= 0) & ccre.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
ccre = ccre[ok].reset_index(drop=True)
print(f"  after edge filter: {len(ccre):,}")

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(ccre), size=need * 2, replace=False if need * 2 <= len(ccre) else True)
    for i in idx:
        if len(pool) >= N:
            break
        r = ccre.iloc[i]
        seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
        if "N" in seq or len(seq) != L:
            continue
        pool.append(seq)
    attempt += 1
assert len(pool) == N

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
