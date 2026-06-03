"""Exp 018 — DHS short-peaks filter (drop broad peaks).

Hypothesis: sharp, narrow DHS peaks (peak width <300bp) represent focused,
high-quality regulatory hotspots — single TF clusters or single regulatory
elements. Broad peaks (>1kb) often represent super-enhancers or chromatin
domains with heterogeneous content; a 200bp window centered on summit
captures only a fraction of the true regulatory unit.

Method: filter to DHSs whose (end - start) <= 300bp (peak width), uniform
sample 50k from this narrower pool (in addition to the 12-topic exclusion).

If eval_01 > 0.6752 → sharp-peak filter improves training signal
If flat/worse → broad peaks are also informative, filter is harmful

This is FILTER (set-membership), not SELECT (rank/weight) — H7-compatible.
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
SEED = 18
MAX_PEAK_WIDTH = 300

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index (with peak width)...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "start", "end", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
df["width"] = df["end"] - df["start"]
print(f"  no-label pool: {len(df):,}; width median={int(df['width'].median())}, "
      f"p90={int(df['width'].quantile(0.9))}")
df = df[df["width"] <= MAX_PEAK_WIDTH].reset_index(drop=True)
print(f"  after width filter (<= {MAX_PEAK_WIDTH}): {len(df):,}")

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
