"""Exp 025 — cCRE class-balanced (25k) + DHS no-label uniform (25k).

Exp 024 (cCRE class-balanced) = 0.6921 NEW BEST. Exp 005 (DHS no-label) = 0.6752.
Prior mixing experiment (exp 006: DHS no-label + random) HURT severely because
random is pure noise. This is a different mix: two HIGH-QUALITY sources.

Hypothesis: DHS no-label adds biosample-derived chromatin context that cCRE's
catalog (which aggregates across all biosamples into a single class label)
doesn't capture. If two good sources span more regulatory grammar than one,
we should beat 0.6921.

Risk: if cCRE class-balanced already covers the space that DHS adds, mixing
just dilutes — flat or hurt. Prior experiments show that single source >
mixed source when one is strictly worse.

Implementation: 25k cCRE class-balanced (~3125 per class) + 25k DHS no-label
uniform. De-duplicate by genomic coord (chrom, mid) before sampling.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"

N, L = 50_000, 200
N_CCRE = N // 2
N_DHS = N - N_CCRE
HALF = L // 2
SEED = 25

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chroms_fa = set(fa.keys())
chrom_len = {}

# ---------------- cCRE class-balanced (25k, ~3125 per class) ----------------
print("loading cCREs...")
ccre = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                   names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
ccre["mid"] = (ccre["start"] + ccre["end"]) // 2
ccre["win_start"] = ccre["mid"] - HALF
ccre["win_end"] = ccre["mid"] + HALF
ccre = ccre[ccre["chrom"].isin(chroms_fa)].reset_index(drop=True)
for c in ccre["chrom"].unique():
    chrom_len[c] = len(fa[c])
ok = (ccre["win_start"] >= 0) & ccre.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
ccre = ccre[ok].reset_index(drop=True)
print(f"  cCREs: {len(ccre):,}")

classes = sorted(ccre["class"].unique())
per_class = N_CCRE // len(classes)
remainder = N_CCRE - per_class * len(classes)

ccre_pool = []
ccre_coords = set()
for ci, cls in enumerate(classes):
    quota = per_class + (1 if ci < remainder else 0)
    sub = ccre[ccre["class"] == cls].reset_index(drop=True)
    print(f"  cCRE {cls}: draw {quota}")
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
            coord = (r.chrom, int(r.mid))
            if coord in ccre_coords:
                continue
            seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
            if "N" in seq or len(seq) != L:
                continue
            collected.append(seq)
            ccre_coords.add(coord)
        attempt += 1
    assert len(collected) == quota, f"{cls}: {len(collected)}/{quota}"
    ccre_pool.extend(collected)
print(f"  cCRE total: {len(ccre_pool)}")

# ---------------- DHS no-label (25k) ----------------
print("loading DHSs (no-label only)...")
dhs = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                  usecols=["seqname", "summit", "component"])
dhs = dhs[~dhs["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
dhs = dhs.rename(columns={"seqname": "chrom"})
dhs["mid"] = dhs["summit"].astype(int)
dhs["win_start"] = dhs["mid"] - HALF
dhs["win_end"] = dhs["mid"] + HALF
dhs = dhs[dhs["chrom"].isin(chroms_fa)].reset_index(drop=True)
for c in dhs["chrom"].unique():
    if c not in chrom_len:
        chrom_len[c] = len(fa[c])
ok = (dhs["win_start"] >= 0) & dhs.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
dhs = dhs[ok].reset_index(drop=True)
print(f"  no-label DHSs: {len(dhs):,}")

dhs_pool, attempt = [], 0
while len(dhs_pool) < N_DHS and attempt < 6:
    need = N_DHS - len(dhs_pool)
    idx = rng.choice(len(dhs), size=need * 2, replace=False if need * 2 <= len(dhs) else True)
    for i in idx:
        if len(dhs_pool) >= N_DHS:
            break
        r = dhs.iloc[i]
        coord = (r.chrom, int(r.mid))
        if coord in ccre_coords:
            continue
        seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
        if "N" in seq or len(seq) != L:
            continue
        dhs_pool.append(seq)
    attempt += 1
print(f"  DHS total: {len(dhs_pool)}")

assert len(ccre_pool) + len(dhs_pool) == N

pool = ccre_pool + dhs_pool
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
