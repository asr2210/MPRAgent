"""005_mpra_dataset — random 50k sample from the published Malinois Table_S2 MPRA.

MAJOR FIND: eval_mapping.json shows that 6/14 eval sets are UKBB/GTEx variant
MPRA evals (eval_02, 03, 05, 06, 11, 12). The Table_S2 dataset contains 798k
MPRA-measured UKBB+GTEx sequences. Direct relevance to half the eval sets.

Hypothesis: Sampling from Table_S2 trains a model that natively understands
the variant MPRA distribution. Should out-perform cCRE/DHS-based libraries on
the UKBB/GTEx evals (eval_02, 03, 05, 06, 11, 12). May lose on synthetic
(eval_08) and DHS (eval_10) — but those were poor for everyone anyway.

Design:
- Random 50k from Table_S2, EXCLUDING chromosomes 7, 13, 19, 21, X
  (eval chromosomes — avoid potential train/eval overlap)
- Eligible pool: ~670k sequences
- Sequence column is already 200bp
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))

# Eval chromosomes per eval_mapping.json: chr7, chr13, chr19, chr21, chrX
EVAL_CHROMS = {'7', '13', '19', '21', 'X'}
MPRA_PATH = '/data/users/arao/.private/mpra_exp/data/Table_S2__MPRA_dataset.txt'
N = 50_000
SEED = 0


def main():
    t0 = time.time()
    print('loading Table_S2 MPRA dataset...', flush=True)
    df = pd.read_csv(MPRA_PATH, sep='\t',
                     usecols=['chr', 'sequence'],
                     dtype={'chr': 'string', 'sequence': 'string'})
    print(f'  {len(df):,} rows in {time.time()-t0:.1f}s', flush=True)

    df['chr'] = df['chr'].astype(str)
    mask = ~df['chr'].isin(EVAL_CHROMS)
    pool = df[mask].reset_index(drop=True)
    print(f'  non-eval-chrom pool: {len(pool):,}', flush=True)

    # Length check
    lengths = pool['sequence'].str.len()
    print(f'  length stats: min={lengths.min()} max={lengths.max()} '
          f'modal={lengths.mode().iloc[0]}', flush=True)
    pool = pool[lengths == 200].reset_index(drop=True)
    print(f'  200bp-only pool: {len(pool):,}', flush=True)

    # Reject any with N
    has_n = pool['sequence'].str.contains('N')
    pool = pool[~has_n].reset_index(drop=True)
    print(f'  no-N pool: {len(pool):,}', flush=True)

    # Dedupe
    pool = pool.drop_duplicates(subset='sequence').reset_index(drop=True)
    print(f'  deduped pool: {len(pool):,}', flush=True)

    # Random sample 50k
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(pool), size=N, replace=False)
    seqs = pool.iloc[idx]['sequence'].tolist()
    seqs = [s.upper() for s in seqs]
    assert len(seqs) == N
    assert all(len(s) == 200 for s in seqs)
    assert all(set(s).issubset({'A', 'C', 'G', 'T'}) for s in seqs[:1000])

    out_path = os.path.join(HERE, 'sequences_0.txt')
    with open(out_path, 'w') as f:
        f.write('\n'.join(seqs) + '\n')
    print(f'wrote {out_path} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
