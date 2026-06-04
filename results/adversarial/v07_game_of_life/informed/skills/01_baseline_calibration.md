# Skill: Baseline calibration for v07 evaluator

## What I observed
The v07 evaluator has a **different performance ceiling** than what
`instructions.md`'s baseline table reports.

- `instructions.md` claims top score around **0.72** (`dhs_topic`)
- `strategies.md` (file `v07.md`) shows top score around **0.40** (`gc_50`)
- My experiment 001 (ENCODE cCRE balanced, 50k) scored **0.3919** — confirming
  the `strategies.md` calibration.

**Trust `strategies.md` (or whichever file is symlinked to the run-specific
strategies file at `runs/v07/.../strategies.md`).** It is calibrated to the
actual evaluator. The `instructions.md` baseline table appears to be from a
different (likely larger / longer-trained) evaluator and is not directly
comparable.

## Real baseline ceiling (v07)
| strategy            | eval_01 |
|---------------------|---------|
| gc_50               | 0.3972  |
| random_uniform      | 0.3951  |
| gc_sweep            | 0.3560  |
| dirichlet_compos.   | 0.3424  |
| at_rich             | 0.2987  |
| gc_rich             | 0.2224  |
| homopolymer_rich    | 0.1748  |
| alternating_ry      | 0.1169  |
| dinuc_repeat        | -0.0306 |

## Key implications
1. **GC ~50% wins.** Uniform composition libraries dominate.
2. **Extreme compositions or wide composition spreads hurt.**
3. **Per-cell-type bias is severe** (K562 > HepG2 > SK-N-SH). Lift the weak
   cell types to lift the mean.
4. **eval_08 is consistently hardest** across all strategies; eval_13 is often
   slightly easier. Other 12 evals cluster tightly in pairs.

## How to know my pipeline is correct
A bare cCRE (or DHS-derived) library at GC ~50% should score around 0.39-0.40,
near `random_uniform`. If you get less than 0.30 from a real genomic
library, something is broken (e.g., bad sequences, wrong length, N bases).
