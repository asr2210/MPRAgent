# Skill: Experiment loop mechanics

## Directory layout
```
libraries/NNN_short_name/
  generate.py        # writes sequences_0.txt (exactly 50000 lines × 200 chars)
  sequences_0.txt
  notes.md           # what tested, what happened, why
  result.json        # written by prepare.py
```

## Running an experiment
```bash
python3 libraries/NNN_name/generate.py
time python3 prepare.py libraries/NNN_name/
```
- prepare.py takes ~2 minutes for 50k sequences (1 seed) on this machine.
- Output is appended to `result.json`. Format: 14 eval blocks
  (`eval_01..eval_14`), each with `mean_r`, `k562_r`, `hepg2_r`, `sknsh_r`.

## After each experiment (per instructions)
1. Write `notes.md` in the experiment dir
2. Append a timestamped entry to `notebook.md`
3. Append a row to `results.tsv`
4. Update / create any skills that should outlive this experiment
5. `git add -A && git commit -m "NNN_name: mean_r=X.XXX" && git push`
   (commit-per-experiment, no batching)

## Sequence constraints (enforced — your generate.py must guarantee)
- Exactly 50000 lines in `sequences_0.txt`
- Each line exactly 200 chars
- Only A, C, G, T (case sensitive — uppercase only is the safe default)

## Results.tsv schema (tab-separated)
`experiment   eval_01..eval_14   time_s   description`

## Python environment notes (this machine)
- `/usr/bin/python3` is Python 3.12.
- Available: numpy, pandas, sklearn, pyfaidx, requests, scipy.
- **NOT available**: torch, pybedtools. So no deep-model active learning
  *locally*; we rely on prepare.py as the only model trainer.
- 20 cores, 3.1TB free under /data/users.

## Genome / reference data location
Store under `data/`. Bulk downloads (hg38 FASTA, DHS Index TSVs, etc.) are fine
but keep them out of git (add to `.gitignore` to keep the repo small).
