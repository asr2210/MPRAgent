# Skill: prepare.py behavior

What I've learned about the black-box `prepare.py` from experiments.

## Interface
```
python prepare.py libraries/NNN_name/sequences.txt
```
Requires `sequences.txt` to be exactly 150,000 lines, each exactly 200
characters from {A, C, G, T}.

Writes:
- `libraries/NNN_name/result.json` — `{eval_01..eval_14: {mean_r, k562_r, hepg2_r, sknsh_r}, time_s}`
- `libraries/NNN_name/model.pt` — trained sequence→activity model (~16 MB)

## Runtime
- ~3400–3500 s (≈57 min) wall-time on the single available GPU (NVIDIA GB10).
  Confirmed on experiment 001 (3433.9 s for 150k random sequences).
- GPU stays around 90+% utilisation during training. Don't queue parallel
  prepare.py runs against the same GPU.
- Plan for ~28 h total compute if running 30 experiments back-to-back.

## What it does (inferred — never read the source)
1. Synthesises an MPRA measurement for each of the 150k sequences in K562,
   HepG2, and SK-N-SH.
2. Trains a sequence→activity model from scratch on (sequence, activity)
   pairs from this library.
3. Evaluates that trained model on 14 anonymous held-out test sets of real
   MPRA-measured sequences and reports Pearson r per cell line and the mean.

## Anchoring results: the random baseline
Pure uniform random 200bp DNA (50% GC, i.i.d.), 150k sequences, seed=0:
- **mean across 14 evals = 0.820**
- K562 = 0.811, HepG2 = 0.820, SKNSH = 0.828
- Range across evals: 0.739 (eval_12, hardest) – 0.908 (eval_08, easiest)

This is the **floor to beat**. Any new library should be compared against
this anchor, not against zero.

## Eval-set behaviour (preliminary, n=1)
Some evals report near-identical numbers on random — possible shared data:
- {02, 06, 14}: ≈ 0.872 each
- {01, 05}:    ≈ 0.776 each
- {03, 10}:    similar (0.856, 0.867)

Hard evals on random: 11 (0.762), 12 (0.739). These will be the most
useful discriminators between library designs — improvements on easy evals
may be saturated.

Easy evals on random: 08 (0.908), 09 (0.889), 14 (0.872).

Until proven otherwise, treat each cluster as one signal when judging
small deltas.

## Operational notes
- After each run: append result row to `results.tsv`, write `notes.md`
  inside the library dir, append an entry to `notebook.md`, then
  `git add -A && git commit -m "NNN_name: mean_r=X.XXX" && git push`.
- Don't read or modify `prepare.py`. It's a black box per instructions.md.
- The output `model.pt` is the trained sequence→activity model. You could
  in principle load it and inspect what it learned, but that's a separate
  diagnostic activity.
