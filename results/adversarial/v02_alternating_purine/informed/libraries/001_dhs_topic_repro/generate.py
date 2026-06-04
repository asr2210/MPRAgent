"""001_dhs_topic_repro — reproduce the dhs_topic baseline at 50k, seed 0.

Anchor against the published baseline mean of eval_01 = 0.7232 (5-seed mean).
A single seed result around 0.71-0.73 confirms my pipeline matches.
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
from data_utils import sample_dhs_sequences  # noqa: E402

DATA = os.path.join(ROOT, 'data')
N = 50_000
SEED = 0


def main():
    t0 = time.time()
    print('loading DHS index...', flush=True)
    dhs = pd.read_csv(os.path.join(DATA, 'DHS_Index.txt.gz'), sep='\t',
                      low_memory=False)
    print(f'  {len(dhs):,} rows in {time.time()-t0:.1f}s', flush=True)

    print('loading NMF mixture...', flush=True)
    with gzip.open(os.path.join(DATA, 'DHS_NMF_Mixture.npy.gz')) as f:
        nmf = np.load(f)  # (16, N_elements)
    assert nmf.shape[1] == len(dhs), (nmf.shape, len(dhs))
    print(f'  shape {nmf.shape} in {time.time()-t0:.1f}s', flush=True)

    # Topic-weighted: per-element weight = sum of NMF loadings across topics.
    # Elements with strong loading on any topic get more weight.
    weights = nmf.sum(axis=0)
    print(f'weights: mean={weights.mean():.4f} max={weights.max():.4f} '
          f'nonzero={(weights>0).sum():,}', flush=True)

    seqs = sample_dhs_sequences(dhs, weights, N, seed=SEED)
    assert len(seqs) == N
    out_path = os.path.join(HERE, 'sequences_0.txt')
    with open(out_path, 'w') as f:
        f.write('\n'.join(seqs) + '\n')
    print(f'wrote {out_path} ({N:,} sequences) in {time.time()-t0:.1f}s',
          flush=True)


if __name__ == '__main__':
    main()
