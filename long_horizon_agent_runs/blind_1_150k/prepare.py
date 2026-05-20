"""
prepare.py — MPRA evaluation harness.

Usage:
    python prepare.py libraries/NNN_name/sequences.txt

Reads exactly 150,000 200bp sequences from sequences.txt, evaluates them,
and writes result.json to the experiment directory.

DO NOT MODIFY THIS FILE.
"""

import sys
import os
import json
import time
import subprocess

_PYTHON = "/data/users/arao/mpra_eval_venv/bin/python"

_RUNNER = """
import sys, json, warnings, os, pickle, base64
warnings.filterwarnings('ignore')

seqs_path = sys.argv[1]
ckpt_path = sys.argv[2]

with open(seqs_path) as f:
    sequences = [l.strip() for l in f if l.strip()]

# Load extended test sets from hidden location
_base = base64.b64decode(b'L2RhdGEvdXNlcnMvYXJhby8ubXByYV9ldmFsX2ludGVybmFs').decode()
_ts_path = os.path.join(_base, 'test_sets_extended.pkl')
with open(_ts_path, 'rb') as f:
    _ts_raw = pickle.load(f)

# Build eval sets dict (same logic as build_eval_sets)
_eval_sets = {}
for name, d in _ts_raw.items():
    labels = d.get('ground_truth', d.get('oracle_labels'))
    _eval_sets[name] = (d['sequences'], labels)
    if 'oracle_labels' in d and 'ground_truth' in d:
        _eval_sets[f'{name}_oracle'] = (d['sequences'], d['oracle_labels'])

# Fixed anonymization mapping (shuffled with seed 2026)
import numpy as np
_names = list(_eval_sets.keys())
_rng = np.random.RandomState(2026)
_shuffled = _names.copy()
_rng.shuffle(_shuffled)
_label_map = {f'eval_{i+1:02d}': name for i, name in enumerate(_shuffled)}
_anon_sets = {label: _eval_sets[name] for label, name in _label_map.items()}

# Label and oracle
from mpra_eval._oracle import label_sequences as _label
from mpra_eval._surrogate import train_and_eval as _train

labels  = _label(sequences, batch_size=512)
metrics = _train(sequences, labels, _anon_sets, seed=0,
                 conv_warm_start=False, save_checkpoint=ckpt_path)

# Build output: per eval set, mean + per-cell-type pearson r
out = {}
for label in _label_map:
    if label in metrics:
        m = metrics[label]
        out[label] = {
            'mean_r': round(m['mean_pearson'], 4),
            'k562_r': round(m['pearson']['K562'], 4),
            'hepg2_r': round(m['pearson']['HepG2'], 4),
            'sknsh_r': round(m['pearson']['SKNSH'], 4),
        }

print(json.dumps(out))
"""


def main():
    if len(sys.argv) != 2:
        print("Usage: python prepare.py libraries/NNN_name/sequences.txt")
        sys.exit(1)

    seqs_path = sys.argv[1]
    if not os.path.isfile(seqs_path):
        print(f"Error: {seqs_path} not found")
        sys.exit(1)

    with open(seqs_path) as f:
        sequences = [line.strip() for line in f if line.strip()]

    if len(sequences) != 150_000:
        print(f"Error: expected 150,000 sequences, got {len(sequences)}")
        sys.exit(1)

    bad = [i for i, s in enumerate(sequences)
           if len(s) != 200 or not all(c in 'ACGT' for c in s)]
    if bad:
        print(f"Error: {len(bad)} sequences have invalid length or characters "
              f"(first bad index: {bad[0]})")
        sys.exit(1)

    print(f"Evaluating {len(sequences):,} sequences against 14 anonymous test sets...")
    t0 = time.perf_counter()

    exp_dir   = os.path.dirname(seqs_path)
    ckpt_path = os.path.join(exp_dir, "model.pt")

    proc = subprocess.run(
        [_PYTHON, "-c", _RUNNER, seqs_path, ckpt_path],
        capture_output=True, text=True
    )

    if proc.returncode != 0:
        print("Evaluation failed:")
        print(proc.stderr[-2000:])
        sys.exit(1)

    elapsed = time.perf_counter() - t0

    scores = json.loads(proc.stdout.strip().splitlines()[-1])
    scores["time_s"] = round(elapsed, 1)

    result_path = os.path.join(exp_dir, "result.json")
    with open(result_path, "w") as f:
        json.dump(scores, f, indent=2)

    print(f"\nResults ({elapsed:.0f}s):")
    for label in sorted(scores):
        if label == "time_s":
            continue
        s = scores[label]
        print(f"  {label}: mean={s['mean_r']:.4f}  "
              f"K562={s['k562_r']:.4f}  HepG2={s['hepg2_r']:.4f}  SKNSH={s['sknsh_r']:.4f}")
    print(f"\nResult written to {result_path}")


if __name__ == "__main__":
    main()
