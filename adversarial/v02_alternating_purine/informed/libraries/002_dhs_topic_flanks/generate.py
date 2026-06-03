"""002_dhs_topic_flanks — DHS-topic positives (25k) + paired genomic flanks (25k).

Hypothesis: blind_claude found that paired flanks help (cCRE+flanks = 0.166).
Does the same recipe applied to DHS positives also help, or is the benefit
cCRE-specific? This decouples "source pool" from "recipe".

Design:
- 25,000 DHS positives: weighted by sum of NMF topic loadings
- 25,000 paired flanks: for each positive, draw an offset in [-3000, -1500] U
  [+1500, +3000] from summit, extract 200bp window. Reject if window overlaps
  any DHS element (so it's truly "non-accessible" background) or contains N.
"""
from __future__ import annotations

import gzip
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
from data_utils import sample_dhs_sequences, _fasta  # noqa: E402

DATA = os.path.join(ROOT, 'data')
N_POS = 25_000
N_FLANK = 25_000
WINDOW = 200
SEED = 0
HALF = WINDOW // 2
FLANK_MIN = 1500
FLANK_MAX = 3000


def build_dhs_interval_lookup(dhs: pd.DataFrame):
    """Return per-chrom sorted (start,end) arrays for fast overlap rejection."""
    by_chrom = {}
    for chrom, g in dhs.groupby('seqname', sort=False):
        starts = g['start'].to_numpy()
        ends = g['end'].to_numpy()
        order = np.argsort(starts)
        by_chrom[chrom] = (starts[order], ends[order])
    return by_chrom


def overlaps_dhs(chrom: str, win_start: int, win_end: int, lookup) -> bool:
    if chrom not in lookup:
        return False
    starts, ends = lookup[chrom]
    # binary search for the first interval whose start >= win_end
    j = np.searchsorted(starts, win_end, side='left')
    # check intervals with index < j whose end > win_start
    # walk back until starts[i] < win_start - margin (intervals can be wide)
    # In DHS, intervals are short (~100-1000bp), so a small walkback is safe.
    i = j - 1
    while i >= 0 and starts[i] > win_start - 5000:
        if ends[i] > win_start and starts[i] < win_end:
            return True
        i -= 1
    return False


def sample_flanks(positives: pd.DataFrame, lookup, fa, seed: int,
                  n_target: int) -> list[str]:
    rng = np.random.default_rng(seed + 1)
    out: list[str] = []
    seen: set[str] = set()
    n_attempt = 0
    n_overlap = n_n = n_dup = n_short = 0
    # Iterate in randomized order, multiple passes if needed.
    idx = rng.permutation(len(positives))
    pass_no = 0
    t0 = time.time()
    while len(out) < n_target and pass_no < 5:
        pass_no += 1
        for i in idx:
            if len(out) >= n_target:
                break
            chrom = positives.iat[int(i), positives.columns.get_loc('seqname')]
            summit = int(positives.iat[int(i), positives.columns.get_loc('summit')])
            # try a few offsets
            for _ in range(3):
                sign = 1 if rng.random() < 0.5 else -1
                off = sign * int(rng.integers(FLANK_MIN, FLANK_MAX + 1))
                center = summit + off
                ws, we = center - HALF, center + HALF
                if ws < 0:
                    n_short += 1
                    continue
                if overlaps_dhs(chrom, ws, we, lookup):
                    n_overlap += 1
                    continue
                try:
                    s = str(fa[chrom][ws:we]).upper()
                except Exception:
                    n_short += 1
                    continue
                if len(s) != WINDOW:
                    n_short += 1
                    continue
                if 'N' in s:
                    n_n += 1
                    continue
                if s in seen:
                    n_dup += 1
                    continue
                seen.add(s)
                out.append(s)
                break
            n_attempt += 1
        # reshuffle for the next pass
        idx = rng.permutation(len(positives))
    print(f'  flanks kept {len(out):,}/{n_target}  attempts={n_attempt:,} '
          f'overlap={n_overlap} N={n_n} dup={n_dup} short={n_short} '
          f'in {time.time()-t0:.1f}s', flush=True)
    if len(out) < n_target:
        raise RuntimeError(f'only {len(out)} flanks accepted')
    return out


def main():
    t0 = time.time()
    print('loading DHS index...', flush=True)
    dhs = pd.read_csv(os.path.join(DATA, 'DHS_Index.txt.gz'), sep='\t',
                      low_memory=False)
    print(f'  {len(dhs):,} rows in {time.time()-t0:.1f}s', flush=True)

    print('loading NMF mixture...', flush=True)
    with gzip.open(os.path.join(DATA, 'DHS_NMF_Mixture.npy.gz')) as f:
        nmf = np.load(f)
    assert nmf.shape[1] == len(dhs)
    weights = nmf.sum(axis=0)

    # Stage 1: sample 25k DHS-topic positives.
    print('sampling 25k DHS-topic positives...', flush=True)
    pos_seqs = sample_dhs_sequences(dhs, weights, N_POS, seed=SEED)

    # To get the positives' metadata we need indices, but sample_dhs_sequences
    # returns just strings. Re-derive: redo the sampling with the same seed
    # and capture indices.
    rng = np.random.default_rng(SEED)
    nz = weights > 0
    pool = np.flatnonzero(nz)
    w = weights[nz]; w = w / w.sum()
    n_draw = min(int(N_POS * 1.3), pool.size)
    drawn_idx = rng.choice(pool, size=n_draw, replace=False, p=w)
    # The data_utils function takes the first N_POS that pass filters, but I
    # don't need exact same indices — I just need 25k DHS rows to pair flanks
    # with. So use drawn_idx[:N_POS] (the order is preserved).
    fa = _fasta()
    pos_meta = dhs.iloc[drawn_idx].reset_index(drop=True)

    # Stage 2: sample flanks paired with positives.
    print('building DHS overlap lookup...', flush=True)
    lookup = build_dhs_interval_lookup(dhs)
    print(f'  lookup built in {time.time()-t0:.1f}s', flush=True)

    print('sampling 25k paired flanks...', flush=True)
    flank_seqs = sample_flanks(pos_meta, lookup, fa, seed=SEED, n_target=N_FLANK)

    all_seqs = pos_seqs + flank_seqs
    assert len(all_seqs) == 50_000
    # Shuffle so positives and flanks aren't blocked
    rng2 = np.random.default_rng(SEED + 42)
    rng2.shuffle(all_seqs)

    out_path = os.path.join(HERE, 'sequences_0.txt')
    with open(out_path, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out_path} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
