"""006_blind013_repro — replicate blind_claude's best recipe.

Blind 013_asym_strat_flank got eval_01 = 0.173. Recipe:
- 15k UNIFORM cCRE centers (preserves natural distribution: ~12k dELS + 2.4k pELS
  + 0.6k PLS + 0.5k CTCF + 0.35k DNH3)
- 5k CTCF-only (BOOST from natural ~0.5k)
- 5k DNase-H3K4me3 (BOOST from natural ~0.35k)
- 25k paired flanks (±1500-3000bp from each positive center, no cCRE overlap)

Goal: confirm my pipeline reaches ~0.17 on the known-best recipe. If yes,
use this as my anchor and build variations. If no, my pipeline has subtle
issues to track down.
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

N_UNIFORM = 15_000
N_CTCF_BOOST = 5_000
N_DNH3_BOOST = 5_000
N_POS = N_UNIFORM + N_CTCF_BOOST + N_DNH3_BOOST  # 25k
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


def extract_centers(df, indices, fa, seen):
    """Extract 200bp around the center for the given indices. Skip dups/N."""
    out_seqs, out_meta = [], []
    for i in indices:
        row = df.iloc[int(i)]
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
    return out_meta, out_seqs


def sample_positives(df, seed):
    rng = np.random.default_rng(seed)
    fa = _fasta()
    seen = set()

    # Component 1: 15k uniform across ALL cCREs
    print('  15k uniform cCREs...', flush=True)
    pool = df.reset_index(drop=True)
    oversample = int(N_UNIFORM * 1.3)
    idx = rng.choice(len(pool), size=oversample, replace=False)
    meta_u, seqs_u = extract_centers(pool, idx[:oversample], fa, seen)
    meta_u = meta_u[:N_UNIFORM]
    seqs_u = seqs_u[:N_UNIFORM]
    assert len(seqs_u) == N_UNIFORM, f'uniform: only {len(seqs_u)}'

    # Component 2: 5k CTCF-only boost (sample from CTCF pool, may overlap with uniform)
    print('  5k CTCF-only boost...', flush=True)
    ctcf_pool = df[df['main_class'] == 'CTCF-only'].reset_index(drop=True)
    oversample = int(N_CTCF_BOOST * 1.3)
    idx = rng.choice(len(ctcf_pool), size=min(oversample, len(ctcf_pool)),
                     replace=False)
    meta_c, seqs_c = extract_centers(ctcf_pool, idx, fa, seen)
    meta_c = meta_c[:N_CTCF_BOOST]
    seqs_c = seqs_c[:N_CTCF_BOOST]
    assert len(seqs_c) == N_CTCF_BOOST, f'CTCF: only {len(seqs_c)}'

    # Component 3: 5k DNase-H3K4me3 boost
    print('  5k DNase-H3K4me3 boost...', flush=True)
    dnh_pool = df[df['main_class'] == 'DNase-H3K4me3'].reset_index(drop=True)
    oversample = int(N_DNH3_BOOST * 1.3)
    idx = rng.choice(len(dnh_pool), size=min(oversample, len(dnh_pool)),
                     replace=False)
    meta_d, seqs_d = extract_centers(dnh_pool, idx, fa, seen)
    meta_d = meta_d[:N_DNH3_BOOST]
    seqs_d = seqs_d[:N_DNH3_BOOST]
    assert len(seqs_d) == N_DNH3_BOOST, f'DNH3: only {len(seqs_d)}'

    meta = meta_u + meta_c + meta_d
    seqs = seqs_u + seqs_c + seqs_d
    print(f'  total positives: {len(seqs):,}', flush=True)
    return meta, seqs


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
        raise RuntimeError(f'only {len(out)}/{n_target} flanks')
    return out


def main():
    t0 = time.time()
    print('loading cCREs...', flush=True)
    ccres = load_ccres()
    print(f'  {len(ccres):,} cCREs', flush=True)

    print('sampling 25k positives (blind 013 recipe)...', flush=True)
    pos_meta, pos_seqs = sample_positives(ccres, SEED)

    print('building cCRE overlap lookup...', flush=True)
    lookup = build_overlap_lookup(ccres)

    print('sampling 25k paired flanks...', flush=True)
    fa = _fasta()
    flank_seqs = sample_flanks(pos_meta, lookup, fa, SEED, N_FLANK)

    all_seqs = pos_seqs + flank_seqs
    assert len(all_seqs) == 50_000
    np.random.default_rng(SEED + 99).shuffle(all_seqs)

    out = os.path.join(HERE, 'sequences_0.txt')
    with open(out, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
