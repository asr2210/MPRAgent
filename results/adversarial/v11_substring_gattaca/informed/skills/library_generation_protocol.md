# Library generation protocol (best practice)

## Constraints
- Exactly 50,000 lines in `sequences_0.txt`
- Each line exactly 200 chars from {A, C, G, T}
- No trailing whitespace issues
- prepare.py is a black box — never modify or read it

## Validated optimal recipe (mean_r = 0.8627, eval_01 = 0.8808)
1. Generate IID Uniform random with seed=2 (best of 3 tested seeds)
2. Rejection-sample to GC ∈ [90, 110] (~67% acceptance)
3. Per-position A↔T and C↔G swap balance
4. Shuffle
5. Write line by line

See `gc_band_rejection.md` and `per_position_balance.md` for code.

## What does NOT help
- Biology (chr22, MPRA, DHS): -0.04 to -0.13 (eval punishes biology)
- Planted TF motifs: -0.027
- Markov dinucleotide: -0.063
- K-mer entropy selection (high or low): ±0.011
- RC pair augmentation: 0
- Per-base count filter: 0
- Homopolymer rejection: 0 (filter too rare to act)
- Gaussian GC weighting: -0.015 (too peaked)
- Uniform GC on band: -0.005 (loses Binomial shape)
- Tighter GC bands [99,101] or [93,107]: -0.05 to -0.01

## Practical tips
- Use `np.random.default_rng(SEED)` for reproducibility
- Batch generation 100k at a time — np vectorization is fast
- Verify per-seq GC range and per-position counts before writing
- Always shuffle before writing — prepare.py may not internally shuffle
- Library generation takes ~30s. prepare.py eval takes ~50s for 50k.
- Budget: ~80s wall per experiment, so 30 experiments = 40 min compute.

## File layout per experiment
```
libraries/NNN_name/
├── generate.py
├── sequences_0.txt
├── result.json     # written by prepare.py
└── notes.md
```

## Commit after each experiment
```bash
git add -A && git commit -m "NNN_name: mean_r=X.XXX"
```
