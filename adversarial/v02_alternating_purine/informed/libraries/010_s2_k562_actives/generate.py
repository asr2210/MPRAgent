"""010_s2_k562_actives — anchor library with K562-active Table_S2 sequences.

Theory v7: Each source contributes signal AND drag per cell-type head.
CTCF dragged K562 because CTCF cCREs have low K562 activity. To raise K562
correlation we need sequences with HIGH K562 activity that span the activity
dynamic range. Table_S2 has measured K562_log2FC -> stratify directly.

Composition (50k) — extends exp 007's winning recipe:
- 10k cCRE uniform (exp 007 baseline)
- 4k cCRE CTCF-only
- 4k cCRE DNase-H3K4me3
- 3k Table_S2 K562-active (top quartile by K562_log2FC, non-eval chroms) — NEW
- 2k Table_S2 random non-eval (was 5k)
- 4k DHS-topic
- 3k synth
- 20k paired flanks
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
from data_utils import _fasta  # noqa: E402

DATA = os.path.join(ROOT, 'data')
WINDOW = 200
HALF = WINDOW // 2
SEED = 0
FLANK_MIN = 1500
FLANK_MAX = 3000
EVAL_CHROMS_BARE = {'7', '13', '19', '21', 'X'}
MPRA_PATH = '/data/users/arao/.private/mpra_exp/data/Table_S2__MPRA_dataset.txt'

N_CCRE_UNI = 10_000
N_CTCF = 4_000
N_DNH3 = 4_000
N_S2_K562 = 3_000  # NEW: K562-top-quartile
N_S2_RAND = 2_000
N_DHS = 4_000
N_SYNTH = 3_000
N_FLANK = 20_000
assert (N_CCRE_UNI + N_CTCF + N_DNH3 + N_S2_K562 + N_S2_RAND
        + N_DHS + N_SYNTH + N_FLANK) == 50_000


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
    seqs, meta = [], []
    for i in indices:
        row = df.iloc[int(i)]
        ws, we = row['center'] - HALF, row['center'] + HALF
        try:
            s = str(fa[row['chrom']][ws:we]).upper()
        except Exception:
            continue
        if len(s) != WINDOW or 'N' in s or s in seen:
            continue
        seen.add(s); seqs.append(s); meta.append((row['chrom'], int(row['center'])))
    return meta, seqs


def sample_class(df, klass, n, seed, seen):
    pool = df if klass is None else df[df['main_class'] == klass]
    pool = pool.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    fa = _fasta()
    idx = rng.choice(len(pool), size=min(int(n * 1.4), len(pool)), replace=False)
    meta, seqs = extract_centers(pool, idx, fa, seen)
    return meta[:n], seqs[:n]


def _load_s2():
    df = pd.read_csv(MPRA_PATH, sep='\t',
                     usecols=['chr', 'sequence', 'K562_log2FC'],
                     dtype={'chr': 'string', 'sequence': 'string',
                            'K562_log2FC': 'float32'})
    df['chr'] = df['chr'].astype(str)
    df = df[~df['chr'].isin(EVAL_CHROMS_BARE)]
    df = df[df['sequence'].str.len() == 200]
    df = df[~df['sequence'].str.contains('N', regex=False)]
    df = df.dropna(subset=['K562_log2FC'])
    df = df.drop_duplicates('sequence').reset_index(drop=True)
    return df


def sample_s2_k562_top(seed, seen, n):
    df = _load_s2()
    cutoff = df['K562_log2FC'].quantile(0.75)
    pool = df[df['K562_log2FC'] >= cutoff].reset_index(drop=True)
    print(f'  K562-top-quartile pool: {len(pool):,} (>= {cutoff:.3f} log2FC)', flush=True)
    rng = np.random.default_rng(seed + 17)
    idx = rng.choice(len(pool), size=min(int(n * 1.3), len(pool)), replace=False)
    out = []
    for i in idx:
        s = pool.iloc[int(i)]['sequence'].upper()
        if s in seen:
            continue
        seen.add(s); out.append(s)
        if len(out) >= n:
            break
    assert len(out) == n, f'K562-top {len(out)}/{n}'
    return out


def sample_s2_random(seed, seen, n):
    df = _load_s2()
    rng = np.random.default_rng(seed + 3)
    idx = rng.choice(len(df), size=int(n * 1.5), replace=False)
    out = []
    for i in idx:
        s = df.iloc[int(i)]['sequence'].upper()
        if s in seen:
            continue
        seen.add(s); out.append(s)
        if len(out) >= n:
            break
    assert len(out) == n
    return out


def sample_dhs(seed, seen, n):
    dhs = pd.read_csv(os.path.join(DATA, 'DHS_Index.txt.gz'), sep='\t', low_memory=False)
    with gzip.open(os.path.join(DATA, 'DHS_NMF_Mixture.npy.gz')) as f:
        nmf = np.load(f)
    weights = nmf.sum(axis=0)
    rng = np.random.default_rng(seed + 5)
    nz = weights > 0; pool = np.flatnonzero(nz)
    w = weights[nz]; w = w / w.sum()
    fa = _fasta()
    idx = rng.choice(pool, size=min(int(n * 1.5), pool.size), replace=False, p=w)
    seqs, meta = [], []
    s_col = list(dhs.columns).index('seqname')
    sm_col = list(dhs.columns).index('summit')
    for i in idx:
        if len(seqs) >= n:
            break
        ch = dhs.iat[int(i), s_col]; su = int(dhs.iat[int(i), sm_col])
        try:
            s = str(fa[ch][su - HALF:su + HALF]).upper()
        except Exception:
            continue
        if len(s) != WINDOW or 'N' in s or s in seen:
            continue
        seen.add(s); seqs.append(s); meta.append((ch, su))
    assert len(seqs) == n
    return meta, seqs


def sample_synth(seed, seen, n):
    rng = np.random.default_rng(seed + 7)
    bases = np.array(['A', 'C', 'G', 'T'])
    out = []
    while len(out) < n:
        s = ''.join(bases[rng.integers(0, 4, size=WINDOW)])
        if s in seen:
            continue
        seen.add(s); out.append(s)
    return out


def sample_flanks(positives, lookup, fa, seed, n, seen):
    rng = np.random.default_rng(seed + 11)
    out = []
    idx = rng.permutation(len(positives))
    pass_no = 0; t0 = time.time()
    while len(out) < n and pass_no < 8:
        pass_no += 1
        for i in idx:
            if len(out) >= n:
                break
            chrom, center = positives[int(i)]
            for _ in range(4):
                sign = 1 if rng.random() < 0.5 else -1
                off = sign * int(rng.integers(FLANK_MIN, FLANK_MAX + 1))
                c2 = center + off
                ws, we = c2 - HALF, c2 + HALF
                if ws < 0 or overlaps(chrom, ws, we, lookup):
                    continue
                try:
                    s = str(fa[chrom][ws:we]).upper()
                except Exception:
                    continue
                if len(s) != WINDOW or 'N' in s or s in seen:
                    continue
                seen.add(s); out.append(s); break
        idx = rng.permutation(len(positives))
    print(f'  flanks: {len(out):,}/{n} ({time.time()-t0:.1f}s)', flush=True)
    if len(out) < n:
        raise RuntimeError(f'only {len(out)}/{n} flanks')
    return out


def main():
    t0 = time.time()
    seen = set()
    print('loading cCREs...', flush=True)
    ccres = load_ccres()

    print(f'cCRE uniform {N_CCRE_UNI}...', flush=True)
    m_u, s_u = sample_class(ccres, None, N_CCRE_UNI, SEED, seen)
    print(f'cCRE CTCF {N_CTCF}...', flush=True)
    m_c, s_c = sample_class(ccres, 'CTCF-only', N_CTCF, SEED + 1, seen)
    print(f'cCRE DNH3 {N_DNH3}...', flush=True)
    m_d, s_d = sample_class(ccres, 'DNase-H3K4me3', N_DNH3, SEED + 2, seen)
    print(f'Table_S2 K562-top {N_S2_K562}...', flush=True)
    s_s2k = sample_s2_k562_top(SEED, seen, N_S2_K562)
    print(f'Table_S2 random {N_S2_RAND}...', flush=True)
    s_s2r = sample_s2_random(SEED, seen, N_S2_RAND)
    print(f'DHS-topic {N_DHS}...', flush=True)
    m_dh, s_dh = sample_dhs(SEED, seen, N_DHS)
    print(f'synth {N_SYNTH}...', flush=True)
    s_sy = sample_synth(SEED, seen, N_SYNTH)

    lookup = build_overlap_lookup(ccres)
    fa = _fasta()
    pos_for_flanks = m_u + m_c + m_d + m_dh
    print(f'flanks {N_FLANK} from {len(pos_for_flanks)} positives...', flush=True)
    s_f = sample_flanks(pos_for_flanks, lookup, fa, SEED, N_FLANK, seen)

    all_seqs = s_u + s_c + s_d + s_s2k + s_s2r + s_dh + s_sy + s_f
    assert len(all_seqs) == 50_000, len(all_seqs)
    np.random.default_rng(SEED + 99).shuffle(all_seqs)

    out = os.path.join(HERE, 'sequences_0.txt')
    with open(out, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
