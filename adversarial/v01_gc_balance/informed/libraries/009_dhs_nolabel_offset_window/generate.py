"""Exp 009 — DHS from non-label-aligned topics with random offset window.

Tests positional augmentation: instead of summit ± 100 (regulatory element always
at position 100), each DHS contributes a window with random offset of -50..+50
from summit. The element appears anywhere from position 50 to 150 within the
200bp window. Forces position-invariant learning.

If this lifts above 0.6752 → positional rigidity was a real limitation.
If flat → model is already position-invariant; need different angle.
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
OFFSET = 50  # ± window offset from summit
SEED = 9

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

# Pre-filter: summit ± (HALF + OFFSET) must fit
df["min_start"] = df["summit"] - HALF - OFFSET
df["max_end"] = df["summit"] + HALF + OFFSET
ok = (df["min_start"] >= 0) & df.apply(lambda r: r["max_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  pool after edge filter: {len(df):,}")

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need = N - len(pool)
    idx = rng.choice(len(df), size=need * 2, replace=False if need * 2 <= len(df) else True)
    offsets = rng.integers(-OFFSET, OFFSET + 1, size=len(idx))
    for k, i in enumerate(idx):
        if len(pool) >= N:
            break
        r = df.iloc[i]
        center = r.summit + offsets[k]
        s, e = center - HALF, center + HALF
        seq = str(fa[r.seqname][s:e]).upper()
        if "N" in seq or len(seq) != L:
            continue
        pool.append(seq)
    attempt += 1
assert len(pool) == N

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
