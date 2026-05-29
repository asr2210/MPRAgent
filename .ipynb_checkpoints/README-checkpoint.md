# MPRAgent: Autonomous AI-Guided MPRA Library Design

## Overview

Massively parallel reporter assays (MPRAs) measure how hundreds of thousands of DNA sequences drive gene expression simultaneously, but their utility for training sequence-to-function models depends critically on which sequences are included. We asked whether an autonomous LLM-based agent could discover, without human guidance, the library design principles that determine how well a surrogate model trained on that library generalises to held-out genomic sequences.

We deployed Claude Opus 4.7 as an autonomous research agent in four independent runs under two conditions:

| Condition | Runs | Description |
|-----------|------|-------------|
| **Blind** | Blind 1, Blind 2 | Agent has no information about existing baselines or strategies |
| **Informed** | Informed 1, Informed 2 | Agent can read a summary of human-designed baseline strategies |

Each run conducted 30 sequential design-evaluate-reflect cycles. In each cycle, the agent:
1. Read its lab notebook and prior results
2. Searched the literature for relevant evidence
3. Designed a 50,000-sequence, 200bp DNA library and wrote a `generate.py`
4. Evaluated the library using a black-box MPRA simulator (`prepare.py`) that returns Pearson r vs ground-truth MPRA measurements across 14 anonymous eval sets
5. Updated its notebook with results and revised its running theory

The primary evaluation metric is **eval_01**: Pearson r of surrogate model predictions vs empirical MPRA ground-truth measurements on a held-out genomic test set (chr7 + chr13, 60,055 sequences). 13 additional evals are included to span different regions of sequence space.

### Key results

| Design | eval_01 (Pearson r vs ground truth) |
|--------|--------------------------------------|
| Random floor (synthetic sequences) | 0.684 |
| Human best (DHS topic-weighted, 50k) | 0.723 |
| Agent best — Blind 1 (026_pels_h3k4me3_combo) | 0.738 |
| Agent best — Blind 2 (013_cCRE_extreme_upweight) | 0.748 |
| Agent best — Informed 1 (010_ccre_iid_human_chicken) | **0.760** |
| Agent best — Informed 2 (015_dhs_70_30_gc_stratified) | 0.751 |

---

## Eval set identities

The 14 eval sets were anonymized during agent runs to prevent gaming. Their identities:

| ID | Named set | Description | Sequences | Labels |
|----|-----------|-------------|-----------|--------|
| **eval_01** | chr7_13_gt | chr7 + chr13 genomic — **primary metric** | 60,055 | Ground-truth MPRA log2FC |
| eval_02 | ukbb_gtex_gt_oracle | UKBB+GTEx variant regions | 59,084 | Oracle labels |
| eval_03 | ukbb_gtex_one_oracle | UKBB+GTEx single-variant regions | 30,505 | Oracle labels |
| eval_04 | chr19_21_X_gt | chr19 + chr21 + chrX genomic | 56,340 | Ground-truth MPRA log2FC |
| eval_05 | ukbb_gtex_gt | UKBB+GTEx variant regions | 59,084 | Ground-truth log2FC |
| eval_06 | ukbb_gtex_both_oracle | UKBB+GTEx both-allele regions | 62,966 | Oracle labels |
| eval_07 | sei_chr7_13 | SEI class regions on chr7+13 | 20,000 | Oracle labels |
| eval_08 | synthetic | Uniform random sequences | 20,000 | Oracle labels |
| eval_09 | chr19_21_X_gt_oracle | chr19 + chr21 + chrX genomic | 56,340 | Oracle labels |
| eval_10 | dhs_chr7_13 | DHS elements on chr7+13 | 20,000 | Oracle labels |
| eval_11 | ukbb_gtex_both | UKBB+GTEx both-allele regions | 62,966 | Ground-truth log2FC |
| eval_12 | ukbb_gtex_one | UKBB+GTEx single-variant regions | 30,505 | Ground-truth log2FC |
| eval_13 | genomic_chr7_13 | Random genomic windows on chr7+13 | 20,000 | Oracle labels |
| eval_14 | chr7_13_gt_oracle | chr7 + chr13 genomic | 60,055 | Oracle labels |

eval_01 is the only eval set combining empirical MPRA ground-truth labels with a genomic test set large enough to be reliable, from held-out chromosomes (chr7+13). eval_04 provides an independent cross-chromosome ground-truth check. eval_05/11/12 assess generalisation to variant-effect prediction.

---

## Baseline strategies

The 14 human-designed baselines span a range of compositional strategies, each evaluated at 7 library sizes (10k–300k sequences) across 5 random seeds:

| Strategy | Description |
|----------|-------------|
| `synth_oracle` | Uniform random sequences (random floor) |
| `dhs_random` | Random DHS elements (Meuleman 2020, 3M elements) |
| `dhs_topic` | DHS elements weighted by NMF topic loadings (**human best at 50k**) |
| `dhs_stratified` | DHS elements stratified by NMF component |
| `sei_random` | Random SEI class regions (Chen 2022, 3M elements) |
| `sei_class` | SEI elements class-balanced across 40 functional classes |
| `dhs_sei` | 50% DHS + 50% SEI |
| `dhs_synth` | 50% DHS + 50% random synthetic |
| `sei_synth` | 50% SEI + 50% random synthetic |
| `dhs_sei_synth` | Equal thirds: DHS + SEI + synthetic |
| `dhs_stratified_sei` | DHS stratified + SEI |
| `dhs_stratified_sei_synth` | DHS stratified + equal SEI + synthetic |
| `mpra_oracle` | Malinois training sequences with oracle labels (ceiling) |
| `mpra_real` | Malinois training sequences with real MPRA labels |

---

## Repository structure

```
MPRAgent/
├── README.md                  # This file
├── instructions.md            # Agent task prompt (given to Claude at start of each run)
├── prepare.py                 # Black-box evaluation harness (read-only for agents)
├── setup.sh                   # Clones boda2, downloads Malinois weights, builds eval pkl
│
├── boda2/                     # Git submodule: CODA/Malinois architecture
│
├── eval/                      # Evaluation infrastructure
│   ├── oracle.py              # Malinois inference oracle: label_sequences(seqs) → (N,3)
│   └── surrogate.py           # Surrogate model trainer: train_and_eval(seqs, labels, ...)
│
├── data/
│   ├── malinois/              # Pretrained Malinois weights (downloaded by setup.sh)
│   │   └── artifacts/
│   │       └── torch_checkpoint.pt
│   └── eval_sets/             # Held-out evaluation sequences and labels (plain text)
│       ├── chr7_13_gt_sequences.txt       # 60,055 sequences (eval_01, primary metric)
│       ├── chr7_13_gt_labels_gt.tsv       # Ground-truth MPRA log2FC [K562, HepG2, SKNSH]
│       └── ...                            # 8 additional named sets (gt + oracle variants)
│
├── agent_runs/                # All four agent run records
│   ├── blind_1/               # Branch 42226_1 — Blind condition, run 1
│   │   ├── results.tsv        # One row per experiment: eval_01..eval_14 + description
│   │   ├── notebook.md        # Full append-only lab notebook
│   │   ├── skills/            # Reusable skills the agent wrote during the run
│   │   └── libraries/
│   │       ├── 001_uniform_random/
│   │       │   ├── generate.py     # Code that built this library
│   │       │   ├── result.json     # prepare.py output (14 eval scores)
│   │       │   └── notes.md        # Agent's post-experiment notes
│   │       └── ...                 # 30 experiments total
│   ├── blind_2/               # Branch 42526_1
│   ├── informed_1/            # Branch 42326_1
│   └── informed_2/            # Branch 42326_2
│
├── baselines/                 # Human-designed baseline strategies (14 strategies × 7 sizes × 5 seeds)
│   ├── run_baselines.py       # Dispatch script
│   ├── run_one.py             # Single-seed worker
│   └── results/
│       └── baselines.csv      # Full results table
│
└── analysis/                  # Code to reproduce all paper figures
    ├── fig2_and_S1_baselines.ipynb    # Fig 2 (composition strip+heatmap) + Fig S1 (learning curves)
    ├── fig3_and_fig4_agent.ipynb      # Fig 3 (agent strip+split) + Fig 4 (generalization)
    └── data/                          # Pre-computed results for figure generation
        ├── baselines_results.csv
        ├── strategies.md
        ├── 42226_1_blind_results.tsv
        ├── 42526_1_blind_results.tsv
        ├── 42326_1_informed_results.tsv
        └── 42326_2_informed_results.tsv
```

---
## The oracle and surrogate

**Oracle** (`eval/oracle.py`): Wraps the pretrained Malinois model, a CNN with 4.1M parameters. Input: list of 200bp DNA strings. Output: (N, 3) float32 log2FC predictions [K562, HepG2, SK-N-SH]. 

```python
from eval.oracle import label_sequences
preds = label_sequences(sequences)  # preds.shape == (N, 3)
```

**Surrogate** (`eval/surrogate.py`): Trains a fresh CNN with the same architecture as Malinois from scratch on oracle-labeled sequences. Same architecture as Malinois; locked training config (Adam, CosineAnnealingWarmRestarts, early stopping at patience=20). ~2 minutes per 50k sequences on a GB10. Bit-exact reproducible given the same seed.

```python
from eval.surrogate import train_and_eval
metrics = train_and_eval(sequences, labels, test_sets={'chr7_13': (seqs, labels)}, seed=0)
```

The surrogate trains from scratch so that generalixation must come entirely from the training sequences, not the architecture.

## Running your own agent (sandbox)

To run your own agent in the same experimental framework used in the paper, copy the sandbox components into a fresh git repository and point an LLM agent at `instructions.md`. Prompt with "Read instructions.md and get started."

### 1. Create a sandbox repo

```bash
mkdir my_mpra_agent_run
cd my_mpra_agent_run
git init

cp -r /path/to/MPRAgent/eval      ./
cp -r /path/to/MPRAgent/data      ./
cp    /path/to/MPRAgent/prepare.py ./
cp    /path/to/MPRAgent/setup.sh   ./
cp    /path/to/MPRAgent/instructions.md ./
```

### 2. Run setup

```bash
bash setup.sh
```

`setup.sh` does four things:
1. Clones the `boda2` repository into `boda2/`
2. Downloads the pretrained Malinois weights (~700 MB) from the public GCS bucket
3. Converts the plain-text eval set files in `data/eval_sets/` to `data/eval_sets.pkl`
4. Runs a sanity check

Requires: Python 3.12, PyTorch 2.7, and `curl`, `wget`, or `gsutil` for the checkpoint download.

```bash
python3.12 -m venv venv
source venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cu128
pip install numpy scipy pandas matplotlib lightning
bash setup.sh
```

### 3. Point your agent at the task

Tell your agent: **"Read `instructions.md` and get started."**

The agent will:
- Create `libraries/NNN_description/` for each experiment
- Write `generate.py` (produces 50,000 × 200bp sequences with seeds 0, 1, 2)
- Run `python prepare.py libraries/NNN_description/` to evaluate
- Append results to `results.tsv` and `notebook.md`
- Commit after each experiment

### What `prepare.py` does

`prepare.py` is a black box from the agent's perspective-- essentially a wet-lab collaborator. For each of the three seed files:
1. Labels 50,000 sequences using the Malinois oracle (`eval/oracle.py`)
2. Trains a surrogate model from scratch on the labeled sequences (`eval/surrogate.py`)
3. Evaluates the surrogate on 14 held-out eval sets and returns Pearson r

Scores are averaged across seeds and written to `libraries/NNN_description/result.json`.

**The agent is not told what the eval sets contain.** The names `eval_01` through `eval_14` are fixed but their identities are not revealed, forcing the agent to discover what makes a good library through experimentation.

Each full evaluation (3 seeds × 50k sequences) takes approximately 15–20 minutes on a modern GPU.

---

## Setup for the full repo (figures + baselines)

To reproduce analyses or re-run baselines, install from the repo root with the boda2 submodule:

```bash
git clone --recurse-submodules git@github.com:asr2210/MPRAgent.git
cd MPRAgent

python3.12 -m venv venv
source venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cu128
pip install numpy scipy pandas matplotlib lightning
bash setup.sh
```

All figures are generated from pre-computed results in `analysis/data/` — no GPU or re-training needed:


## Citation

The Malinois model used as oracle (cite as described in the boda2 repository):

```bibtex
@software{sjgosai_boda2,
  author  = {Gosai, Sager J},
  title   = {boda2: Computational Optimization of DNA Activity},
  url     = {https://github.com/sjgosai/boda2},
  year    = {2024},
}
```

---

## License

Code in this repository is released under the MIT License. The Malinois model weights are subject to the boda2 license terms. See `boda2/LICENSE` and `boda2/LICENSE.mit`.
