"""
Experiment 012 — ccre_heavy_mix

008 (50/50 DHS+cCRE) is current best at 0.5671.
011 (cell-type-targeted DHS + cCRE) is 0.5688 — tie.

This experiment sweeps the DHS:cCRE ratio toward cCRE-heavy:
  15 k DHS uniform + 35 k cCRE class-balanced.

Logic:
- cCRE class-balanced alone (006) = 0.5637; DHS uniform alone (002) =
  0.5627. Both very similar standalone.
- Mix (008) gains +0.003 over either.
- If the gain is from cCRE's class-balanced diversity, more cCRE =
  more gain. Diminishing-returns inflection is the interesting test.
- Class-balanced cCRE forces over-representation of rare classes
  (CA-CTCF, TF-only); doubling cCRE roughly doubles their absolute
  count, which could expose more rare-class regulatory grammar.

If 012 > 011: cCRE-heavy direction is real, push further.
If 012 < 008: 50/50 is closer to optimal; cCRE has saturating return.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SEED = 0
N = 50_000
N_DHS = 15_000
N_CCRE = 35_000
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def main():
    rng = np.random.default_rng(SEED)
    # DHS uniform
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    dhs_idx = rng.choice(dhs_seqs.shape[0], size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    # cCRE class-balanced
    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    classes = sorted(ccre_meta["cls"].unique())
    per_cls = N_CCRE // len(classes)  # 35000/8 = 4375
    parts = []
    for c in classes:
        pool = np.where(ccre_meta["cls"].to_numpy() == c)[0]
        n = min(per_cls, len(pool))
        sel = rng.choice(pool, size=n, replace=False) if len(pool) >= n else pool.copy()
        parts.append(sel)
        print(f"  cCRE {c}: pool={len(pool):,}, sampled={len(sel):,}")
    ccre_idx = np.concatenate(parts)
    if len(ccre_idx) < N_CCRE:
        used = set(ccre_idx.tolist())
        extra_pool = np.array([i for i in np.where(ccre_meta["cls"].to_numpy() == "dELS")[0]
                               if i not in used])
        ccre_idx = np.concatenate(
            [ccre_idx, rng.choice(extra_pool, size=N_CCRE - len(ccre_idx), replace=False)]
        )
    ccre_rows = np.asarray(ccre_seqs[ccre_idx])

    all_rows = np.concatenate([dhs_rows, ccre_rows], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} ({N_DHS} DHS + {N_CCRE} cCRE)")


if __name__ == "__main__":
    main()
