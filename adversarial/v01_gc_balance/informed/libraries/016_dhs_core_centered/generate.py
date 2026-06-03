"""Exp 016 — DHS centered on `core_midpoint` instead of `summit`.

The DHS Index has `summit` (single bp of max signal) and `core_start`/`core_end`
(consensus core region, a sub-interval of the full peak). For most DHSs the
core is wider than the summit and represents the region consistently
accessible across biosamples. Centering on core_midpoint may anchor the
window more on the "consensus regulatory element" than on a sharp summit
which can be sample-specific.

Method: same filter as exp 005 (12 non-label-aligned topics, uniform sample).
Center each 200bp window on (core_start + core_end) // 2 instead of summit.

If eval_01 > 0.6752 → core-centered is a better anchor.
If flat → summit and core_midpoint anchor similarly enough.
If worse → summit was actually closer to the regulatory hotspot.
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
SEED = 16

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index (with core)...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "core_start", "core_end", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
df = df.dropna(subset=["core_start", "core_end"]).reset_index(drop=True)

df["core_mid"] = ((df["core_start"] + df["core_end"]) // 2).astype(int)
shift = (df["core_mid"] - df["summit"]).abs()
print(f"  pool: {len(df):,}; |core_mid - summit| median={int(shift.median())}, "
      f"p90={int(shift.quantile(0.9))}, max={int(shift.max())}")

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = (df["core_mid"] - HALF).astype(int)
df["win_end"] = (df["core_mid"] + HALF).astype(int)
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  after edge filter: {len(df):,}")

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
