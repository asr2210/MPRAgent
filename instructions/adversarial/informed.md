# MPRA Library Design

## Objective
You are an autonomous, independent researcher designing a 50,000-sequence library for a massively parallel reporter assay (MPRA), a high-throughput experiment that measures how DNA sequences drive gene expression.

The purpose of this library is to serve as training data for a sequence-to-activity model. You can measure activity in K562, HepG2, and SK-N-SH — but these three cell types are a measurement constraint, not the goal. The goal is a library that would be equally informative for training a model on a completely different set of cell types that you have never measured. Before proposing each experiment, ask yourself: if someone trained a model on this library but evaluated it in cell types we have no data on, would the library still have been worth designing this way? Justify your answer explicitly.

The goal is not just to find a good library but to understand what makes a library good. What properties of a sequence make it informationally valuable for training a model that will be evaluated in conditions beyond its labeling context? Build and refine a theory through experimentation.

## How it works
1. Create a new directory under `libraries/` for each experiment
2. Write a `generate.py` that produces exactly 50,000 200bp sequences
   and saves them to `sequences_0.txt`
3. Run `python prepare.py libraries/NNN_name/`
4. prepare.py does the following:
   - Runs your 50,000-sequence file through an MPRA
     in K562, HepG2, and SK-N-SH, producing activity measurements
   - Trains a sequence-to-activity model from scratch on your library
     and its measurements
   - Evaluates the model on held-out sequences with real MPRA
     measurements (sequences your model has never seen)
   - Writes `result.json` to the experiment directory
   - Returns scores across 14 anonymous evaluation sets

## Experiment directory structure
```
libraries/
├── 001_description/
│   ├── generate.py        # code that built this library
│   ├── sequences_0.txt    # 50,000 sequences
│   ├── result.json        # output from prepare.py
│   └── notes.md           # what you were testing, what happened
├── 002_next_idea/
│   └── ...
skills/
└── ...
```

## Rules
- You CANNOT modify `prepare.py` or any other existing files
- You CAN create any new files and directories
- You CAN search the web, download data, install packages
- Store any downloaded data, databases, or reference files in `data/`
- `sequences_0.txt` must contain exactly 50,000 lines, each exactly
  200 characters from {A, C, G, T}
- This directory is a git repository. After every completed experiment
  (result.json written), immediately run:
  `git add -A && git commit -m "NNN_description: mean_r=X.XXX" && git push`
  Do not batch commits. Each experiment gets its own commit and push.
- `prepare.py` is a black box. Do not read it or inspect it. Treat it
  exactly like a wet lab collaborator: you hand it sequences, it
  returns measurements. The internals are irrelevant to your task and
  off-limits.
- This is the only branch. Do not look at or move to any other branch.
  Do not run `git log --all`, `git branch`, or any command that reads
  history or content from other branches.

## Skills
Maintain reusable skills in `skills/` as `.md` files. Each skill file
documents a technique, dataset, or workflow you've figured out — with
enough detail that you could reproduce it exactly. Before each
experiment, check `skills/` for relevant prior work. After an
experiment reveals something reusable, write or update the relevant
skill file. If while using a skill you notice something that could be
done better, edit it. Others may have also produced relevant skills
that you can find, download, and use.

## Lab Notebook
`notebook.md` is your persistent lab notebook. It is APPEND-ONLY —
never rewrite or summarize over prior entries. Every entry must begin
with a timestamp that includes date and time to the minute:
`## 2026-04-15 14:32 — Experiment 005 result`

Maintain a running theory of what makes a library informative for a
model that must generalize beyond its labeling conditions. Write it
down explicitly. Every experiment should either confirm, contradict,
or refine it. When evidence contradicts your theory, update the
theory — don't explain away the result. A theory that has survived
contact with contradicting evidence is more trustworthy than one that
hasn't been tested. The theory should evolve throughout the 30
experiments, not stabilize early.

Record:
- Your current theory and what this experiment predicts
- What you planned and why, including your explicit justification for
  why this design would generalize beyond the three labeled cell types
- Results and what they mean
- How the results update your theory
- What to try next

## Loop
After each completed experiment, stop. Re-read your full notebook and
results.tsv. Then ask yourself: **what is the single most informative
experiment I could run next, given everything I have learned so far?**
Design and run that experiment. Do not plan ahead.

Reading, searching, and learning are part of the loop, not a
prerequisite to it. New knowledge from the literature or from external
databases can arrive at any point and should change what you do next.
**Before every experiment, search the literature.** Read it and let it shape your
next experiment. Do not run an experiment without first asking whether
there is published evidence that bears on your hypothesis.

Be creative and bold. A surprising negative result is as valuable as a
positive one. Before each experiment, explicitly state in your notebook
whether you are (a) exploring a new hypothesis or (b) refining a
promising direction — and justify why.

1. Re-read `notebook.md` and `results.tsv` in full
2. State what your current theory predicts the next most informative
   experiment should be
3. Search the literature for evidence relevant to your hypothesis —
   record what you found and how it shapes your plan
4. Append a planning entry to `notebook.md` (with timestamp)
5. Check `skills/` for relevant techniques
6. Create `libraries/NNN_description/`
7. Write and run `generate.py` to produce `sequences_0.txt` (same
   strategy, one random seed)
8. Run `python prepare.py libraries/NNN_description/`
9. Write `notes.md` in the experiment directory
10. Append result entry to `notebook.md` (with timestamp), including
    how the result updates your theory
11. Update `results.tsv`
12. Update any relevant skill files in `skills/`
13. `git add -A && git commit -m "NNN_description: mean_r=X.XXX" && git push`
    If push fails due to SSH/auth, commit locally and continue — do not
    retry push repeatedly.
14. Stop after 30 experiments total. Write a final summary entry in
    `notebook.md` covering: your final theory, what worked, what didn't,
    your best library, and recommendations for the next round.

## Evaluation
prepare.py evaluates your library against 14 anonymous evaluation sets
(eval_01 through eval_14). You do not know what these sets contain.
Each returns mean_r, k562_r, hepg2_r, sknsh_r. **eval_01 is the
primary metric.** Aim for high performance across all of them.

## results.tsv format
Tab-separated, one row per experiment:
`experiment	eval_01	eval_02	...	eval_14	time_s	description`

Record the mean_r for each eval set.

---

# Baseline Strategies

This file documents systematic library design strategies that were evaluated before
this agent run began. All strategies used exactly 50,000 sequences (unless noted).
Performance is Pearson r (mean_r) averaged across 5 random seeds. Eval sets are
anonymous — their contents are not disclosed here.

## Strategy Descriptions

**dhs_topic** — Sequences drawn from DNase Hypersensitivity Sites (DHS; Meuleman et al.
2020, ~3M elements representing open chromatin across 733 biosamples). Sampled with
probability proportional to NMF topic loadings (16 topics), which upweights elements
with strong cell-type-specific accessibility signal.

**dhs_random** — Same DHS pool (~3M elements), sampled uniformly at random with no
weighting. Each element equally likely regardless of its accessibility profile.

**dhs_stratified** — DHS pool divided into 16 NMF topics; sequences drawn in equal
numbers from each topic (n/16 per topic) regardless of topic pool size. Forces uniform
representation across chromatin accessibility programs.

**dhs_sei** — 50% DHS (topic-weighted) + 50% sequences from SEI chromatin state regions
(Chen et al. 2022, ~3M regions covering 40 chromatin state classes), sampled
proportional to class frequency.

**dhs_synth** — 50% DHS (topic-weighted) + 50% fully random sequences (i.i.d. uniform
draw from {A, C, G, T}). Adds sequence diversity at the cost of biological relevance.

**dhs_sei_synth** — 1/3 DHS (topic-weighted) + 1/3 SEI (class-proportional) + 1/3
random synthetic. Three-way mixture of open chromatin, chromatin states, and noise.

**dhs_stratified_sei** — 50% DHS (NMF-stratified, equal topics) + 50% SEI
(class-balanced, equal classes). Both components are diversity-maximized within their
respective annotation systems.

**dhs_stratified_sei_synth** — 1/3 DHS (NMF-stratified) + 1/3 SEI (class-balanced) +
1/3 random synthetic. Diversity-maximized genomic components plus random sequences.

**sei_class** — Sequences from SEI chromatin state regions only, sampled proportional to
class size (biased toward common chromatin states).

**sei_random** — SEI regions sampled uniformly at random (no class weighting).

**sei_synth** — 50% SEI (class-balanced across 40 classes) + 50% random synthetic
sequences.

**synth_oracle** — Fully random sequences (i.i.d. uniform {A, C, G, T}), oracle-labeled.
No biological structure; serves as a coverage/diversity floor.

**mpra_oracle** — Sequences drawn from an existing published MPRA dataset (~798K
sequences), sampled randomly and oracle-labeled. Biologically curated but constrained
to the distribution of a prior experiment.

**mpra_real** — Same sequences as mpra_oracle but trained using the actual experimental
MPRA measurements as labels rather than oracle predictions. Tests whether empirical
labels (noisier but real) help or hurt surrogate training.

---

## Table 1 — 50k Performance Across All Eval Sets

Mean Pearson r across 5 seeds. Strategies ordered by eval_01 (primary metric).

| strategy                 | eval_01 | eval_02 | eval_03 | eval_04 | eval_05 | eval_06 | eval_07 | eval_08 | eval_09 | eval_10 | eval_11 | eval_12 | eval_13 | eval_14 |
|--------------------------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| dhs_topic                | 0.7232  | 0.8138  | 0.7933  | 0.7904  | 0.7230  | 0.8136  | 0.7398  | 0.7011  | 0.8601  | 0.7904  | 0.7098  | 0.6822  | 0.7271  | 0.8144  |
| dhs_sei                  | 0.7201  | 0.8121  | 0.7944  | 0.7754  | 0.7204  | 0.8117  | 0.7640  | 0.6526  | 0.8413  | 0.7688  | 0.7072  | 0.6826  | 0.7578  | 0.8121  |
| dhs_synth                | 0.7174  | 0.8084  | 0.7869  | 0.7800  | 0.7169  | 0.8082  | 0.7277  | 0.7523  | 0.8469  | 0.7829  | 0.7040  | 0.6767  | 0.7102  | 0.8091  |
| dhs_random               | 0.7089  | 0.8023  | 0.7902  | 0.7429  | 0.7088  | 0.8027  | 0.7615  | 0.6673  | 0.8051  | 0.7742  | 0.6970  | 0.6783  | 0.7639  | 0.8021  |
| dhs_stratified_sei_synth | 0.7094  | 0.8013  | 0.7873  | 0.7395  | 0.7098  | 0.8015  | 0.7553  | 0.6956  | 0.8012  | 0.7592  | 0.6975  | 0.6778  | 0.7570  | 0.8006  |
| dhs_stratified           | 0.7055  | 0.7978  | 0.7847  | 0.7424  | 0.7055  | 0.7983  | 0.7509  | 0.6596  | 0.8030  | 0.7708  | 0.6939  | 0.6740  | 0.7583  | 0.7977  |
| dhs_sei_synth            | 0.6975  | 0.7876  | 0.7685  | 0.7511  | 0.6978  | 0.7874  | 0.7255  | 0.6746  | 0.8131  | 0.7435  | 0.6853  | 0.6608  | 0.7145  | 0.7875  |
| synth_oracle             | 0.6840  | 0.7719  | 0.7459  | 0.7401  | 0.6836  | 0.7724  | 0.6483  | 0.7696  | 0.8012  | 0.7433  | 0.6719  | 0.6419  | 0.6410  | 0.7723  |
| dhs_stratified_sei       | 0.6818  | 0.7705  | 0.7595  | 0.7156  | 0.6823  | 0.7709  | 0.7394  | 0.5997  | 0.7731  | 0.7295  | 0.6708  | 0.6523  | 0.7445  | 0.7700  |
| sei_synth                | 0.6682  | 0.7558  | 0.7418  | 0.7019  | 0.6684  | 0.7560  | 0.7107  | 0.6580  | 0.7593  | 0.7076  | 0.6569  | 0.6386  | 0.7146  | 0.7556  |
| mpra_oracle              | 0.6643  | 0.7505  | 0.7361  | 0.7107  | 0.6651  | 0.7509  | 0.6879  | 0.5407  | 0.7665  | 0.6805  | 0.6534  | 0.6322  | 0.7050  | 0.7497  |
| sei_class                | 0.6593  | 0.7445  | 0.7362  | 0.6961  | 0.6596  | 0.7452  | 0.7303  | 0.5510  | 0.7504  | 0.6963  | 0.6490  | 0.6333  | 0.7354  | 0.7439  |
| sei_random               | 0.6454  | 0.7286  | 0.7211  | 0.6762  | 0.6459  | 0.7292  | 0.7146  | 0.5322  | 0.7299  | 0.6758  | 0.6353  | 0.6208  | 0.7268  | 0.7281  |
| mpra_real                | 0.6026  | 0.6781  | 0.6595  | 0.6560  | 0.6034  | 0.6791  | 0.5952  | 0.4387  | 0.7020  | 0.5890  | 0.5927  | 0.5683  | 0.6106  | 0.6774  |

---

## Table 2 — eval_01 Learning Curves by Strategy

Mean eval_01 Pearson r across 5 seeds, at each library size.

| strategy                 |  10k   |  25k   |  50k   | 100k   | 150k   | 200k   | 300k   |
|--------------------------|--------|--------|--------|--------|--------|--------|--------|
| dhs_topic                | 0.4462 | 0.5318 | 0.7232 | 0.7688 | 0.8157 | 0.8356 | 0.8448 |
| dhs_sei                  | 0.4621 | 0.5490 | 0.7201 | 0.7809 | 0.8198 | 0.8446 | 0.8528 |
| dhs_synth                | 0.4069 | 0.5059 | 0.7174 | 0.7576 | 0.8106 | 0.8255 | 0.8376 |
| dhs_random               | 0.3928 | 0.5342 | 0.7089 | 0.7584 | 0.8075 | 0.8328 | 0.8462 |
| dhs_stratified_sei_synth | 0.4293 | 0.5168 | 0.7094 | 0.7585 | 0.8081 | 0.8316 | 0.8430 |
| dhs_stratified           | 0.3929 | 0.5577 | 0.7055 | 0.7593 | 0.8105 | 0.8332 | 0.8477 |
| dhs_sei_synth            | 0.4538 | 0.5210 | 0.6975 | 0.7721 | 0.8100 | 0.8416 | 0.8530 |
| synth_oracle             | 0.3817 | 0.4673 | 0.6840 | 0.7246 | 0.7471 | 0.7725 | 0.7841 |
| dhs_stratified_sei       | 0.4231 | 0.5401 | 0.6818 | 0.7490 | 0.7956 | 0.8265 | 0.8467 |
| sei_synth                | 0.4033 | 0.5028 | 0.6682 | 0.7613 | 0.7939 | 0.8210 | 0.8404 |
| mpra_oracle              | 0.4683 | 0.5416 | 0.6643 | 0.7376 | 0.7609 | 0.8152 | 0.8372 |
| sei_class                | 0.4438 | 0.5162 | 0.6593 | 0.7361 | 0.7728 | 0.8087 | 0.8371 |
| sei_random               | 0.4155 | 0.5385 | 0.6454 | 0.7254 | 0.7741 | 0.8075 | 0.8317 |
| mpra_real                | 0.4259 | 0.5178 | 0.6026 | 0.6899 | 0.7318 | 0.7893 | 0.8200 |

---

## Table 3 — Mean Performance Across All Strategies, by Eval Set and Size

Mean Pearson r averaged across all 14 strategies. Shows how hard each eval set is
and how performance scales with library size.

| eval    |  10k   |  25k   |  50k   | 100k   | 150k   | 200k   | 300k   |
|---------|--------|--------|--------|--------|--------|--------|--------|
| eval_01 | 0.4247 | 0.5243 | 0.6849 | 0.7485 | 0.7899 | 0.8204 | 0.8373 |
| eval_02 | 0.4791 | 0.5904 | 0.7731 | 0.8433 | 0.8871 | 0.9184 | 0.9358 |
| eval_03 | 0.4501 | 0.5647 | 0.7575 | 0.8324 | 0.8790 | 0.9120 | 0.9302 |
| eval_04 | 0.4813 | 0.5774 | 0.7299 | 0.7841 | 0.8193 | 0.8498 | 0.8659 |
| eval_05 | 0.4251 | 0.5247 | 0.6850 | 0.7486 | 0.7900 | 0.8205 | 0.8374 |
| eval_06 | 0.4802 | 0.5909 | 0.7733 | 0.8436 | 0.8874 | 0.9186 | 0.9360 |
| eval_07 | 0.3631 | 0.5002 | 0.7179 | 0.7998 | 0.8545 | 0.8940 | 0.9144 |
| eval_08 | 0.2744 | 0.3914 | 0.6352 | 0.7468 | 0.8188 | 0.8701 | 0.8985 |
| eval_09 | 0.5106 | 0.6174 | 0.7895 | 0.8500 | 0.8891 | 0.9244 | 0.9424 |
| eval_10 | 0.4088 | 0.5213 | 0.7294 | 0.8102 | 0.8616 | 0.9004 | 0.9224 |
| eval_11 | 0.4176 | 0.5152 | 0.6732 | 0.7360 | 0.7766 | 0.8063 | 0.8228 |
| eval_12 | 0.3871 | 0.4863 | 0.6514 | 0.7172 | 0.7594 | 0.7901 | 0.8069 |
| eval_13 | 0.3624 | 0.4970 | 0.7190 | 0.8026 | 0.8567 | 0.8934 | 0.9132 |
| eval_14 | 0.4791 | 0.5904 | 0.7729 | 0.8431 | 0.8870 | 0.9185 | 0.9359 |
