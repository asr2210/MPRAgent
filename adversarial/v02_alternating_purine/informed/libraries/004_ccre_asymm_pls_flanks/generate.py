"""004_ccre_asymm_pls_flanks — asymmetric cCRE + PLS boost + paired flanks.

Hypothesis: Blind 013 (15k uniform + 5k CTCF + 5k DNH3 + 25k flanks) hit
0.173 but starved PLS (only ~570 in 15k uniform). Per LentiMPRA (Jan 2025),
promoters have more shared TF grammar across cell types than enhancers, so
PLS may help cross-cell-type generalization (and eval_13 specifically, which
blind noted dropped in 013).

Design (25k positives + 25k flanks = 50k):
- 10k dELS (dominant enhancer signal; reduced from blind 013's ~12k effective)
- 4k pELS (proximal enhancers — promoter-adjacent grammar)
- 3k PLS (BOOST from natural ~0.5% to 12% — addresses 013's gap)
- 4k CTCF-only (boost from natural ~3%)
- 4k DNase-H3K4me3 (boost from natural ~2%)
- 25k paired flanks (±1500-3000bp from each positive center, no cCRE overlap)
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
from data_utils import _fasta  # noqa: E402

DATA = os.path.join(ROOT, 'data')
WINDOW = 200
HALF = WINDOW // 2
SEED = 0
FLANK_MIN = 1500
FLANK_MAX = 3000

CLASS_TARGETS = {
    'dELS': 10_000,
    'pELS': 4_000,
    'PLS': 3_000,
    'CTCF-only': 4_000,
    'DNase-H3K4me3': 4_000,
}
N_POS = sum(CLASS_TARGETS.values())  # 25_000
N_FLANK = 25_000


def load_ccres():
    df = pd.read_csv(os.path.join(DATA, 'GRCh38-cCREs.bed'), sep='\t',
                     header=None,
                     names=['chrom', 'start', 'end', 'rDHS', 'accession', 'classes'])
    df['main_class'] = df['classes'].str.split(',').str[0]
    df['center'] = ((df['start'] + df['end']) // 2).astype(int)
    return df


def build_overlap_lookup(df):
    out = {}
    for chrom, g in df.groupby('chrom', sort=False):
        starts = g['start'].to_numpy()
        ends = g['end'].to_numpy()
        order = np.argsort(starts)
        out[chrom] = (starts[order], ends[order])
    return out


def overlaps(chrom, ws, we, lookup):
    if chrom not in lookup:
        return False
    starts, ends = lookup[chrom]
    j = np.searchsorted(starts, we, side='left')
    i = j - 1
    while i >= 0 and starts[i] > ws - 5000:
        if ends[i] > ws and starts[i] < we:
            return True
        i -= 1
    return False


def sample_positives(df, targets, seed):
    rng = np.random.default_rng(seed)
    fa = _fasta()
    out_seqs, out_meta = [], []
    seen = set()
    for klass, n_need in targets.items():
        pool = df[df['main_class'] == klass].reset_index(drop=True)
        oversample = min(int(n_need * 1.3), len(pool))
        idx = rng.choice(len(pool), size=oversample, replace=False)
        kept = 0
        for i in idx:
            if kept >= n_need:
                break
            row = pool.iloc[int(i)]
            ws, we = row['center'] - HALF, row['center'] + HALF
            try:
                s = str(fa[row['chrom']][ws:we]).upper()
            except Exception:
                continue
            if len(s) != WINDOW or 'N' in s or s in seen:
                continue
            seen.add(s)
            out_seqs.append(s)
            out_meta.append((row['chrom'], int(row['center'])))
            kept += 1
        print(f'  {klass}: {kept:,}/{n_need:,}', flush=True)
    return out_meta, out_seqs


def sample_flanks(positives, lookup, fa, seed, n_target):
    rng = np.random.default_rng(seed + 11)
    out, seen = [], set()
    n_overlap = n_n = n_dup = n_short = 0
    idx = rng.permutation(len(positives))
    pass_no = 0
    t0 = time.time()
    while len(out) < n_target and pass_no < 6:
        pass_no += 1
        for i in idx:
            if len(out) >= n_target:
                break
            chrom, center = positives[int(i)]
            for _ in range(4):
                sign = 1 if rng.random() < 0.5 else -1
                off = sign * int(rng.integers(FLANK_MIN, FLANK_MAX + 1))
                c2 = center + off
                ws, we = c2 - HALF, c2 + HALF
                if ws < 0:
                    n_short += 1
                    continue
                if overlaps(chrom, ws, we, lookup):
                    n_overlap += 1
                    continue
                try:
                    s = str(fa[chrom][ws:we]).upper()
                except Exception:
                    n_short += 1
                    continue
                if len(s) != WINDOW or 'N' in s or s in seen:
                    if s in seen:
                        n_dup += 1
                    elif 'N' in s:
                        n_n += 1
                    else:
                        n_short += 1
                    continue
                seen.add(s)
                out.append(s)
                break
        idx = rng.permutation(len(positives))
    print(f'  flanks kept {len(out):,}/{n_target} overlap={n_overlap} '
          f'N={n_n} dup={n_dup} short={n_short} in {time.time()-t0:.1f}s',
          flush=True)
    if len(out) < n_target:
        raise RuntimeError(f'only {len(out)}/{n_target} flanks accepted')
    return out


def main():
    t0 = time.time()
    print('loading cCREs...', flush=True)
    ccres = load_ccres()
    print(f'  {len(ccres):,} cCREs in {time.time()-t0:.1f}s', flush=True)

    print(f'sampling {N_POS:,} positives (asymmetric by class)...', flush=True)
    pos_meta, pos_seqs = sample_positives(ccres, CLASS_TARGETS, SEED)
    assert len(pos_seqs) == N_POS

    print('building cCRE overlap lookup...', flush=True)
    lookup = build_overlap_lookup(ccres)

    print(f'sampling {N_FLANK:,} paired flanks...', flush=True)
    fa = _fasta()
    flank_seqs = sample_flanks(pos_meta, lookup, fa, SEED, N_FLANK)

    all_seqs = pos_seqs + flank_seqs
    assert len(all_seqs) == 50_000
    np.random.default_rng(SEED + 99).shuffle(all_seqs)

    out_path = os.path.join(HERE, 'sequences_0.txt')
    with open(out_path, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out_path} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
