# Experiment 015 — Motif-tiled (structured non-random subspace)

## Design
Each 200bp sequence built by tiling JASPAR TF motifs (length 6-15bp) with
1-3bp random spacers, padded to 200bp. Result: high motif density,
fundamentally structured sequence subspace.

GC came out 0.473 (close to 0.5), std=0.044. Library mean composition is
near-50, but the sequence STRUCTURE is fundamentally different from
random or implanted-motif libraries.

## Result
- eval_01 = **0.3770** (Δ vs 007 = -0.020)
- mean_r = **0.3651** (Δ vs 007 = -0.021)
- Per-cell on eval_01: K562=0.592, HepG2=0.418, SK-N-SH=0.122

Drop of 0.018-0.021 is well outside noise floor (0.003). This is a
real, large negative effect.

## Interpretation
**Structured motif-dense subspace is WORSE than random GC-50.** The
0.397 plateau is the OPTIMUM for this evaluator, not a saturation
ceiling. Sequence space outside "GC≈0.5 unstructured ACGT" is penalized.

Why:
- Real MPRA test sequences likely have motif density similar to natural
  enhancers (~3-5 motifs / 200bp), NOT tiled. Training on tiled
  sequences teaches the model to expect dense motif structure that
  doesn't generalize.
- Training-set dinucleotide statistics are far from natural — k-mer
  distribution dominated by motif-internal patterns (e.g., NF-κB consensus
  GGGRNNYYCC creates artificial GC clustering).
- The model probably learns to over-attend to "looks-like-a-motif" rather
  than "looks-like-a-real-enhancer-with-some-motifs."

## Theory v9
The 0.397 plateau is structural: it's the optimal performance of the
(model, training-size=50k, evaluator) pipeline given that training
sequences come from any natural-like distribution. Deviations of the
training subspace AWAY FROM natural-like distribution (motif-tiled,
extreme composition, etc.) cost 0.01-0.02.

Combined with 014's noise floor:
- Within natural-like subspace (GC≈0.5 random or near-random): plateau
  at 0.397 ± 0.003 (training noise)
- Outside natural-like subspace: significant penalty

The "art" of library design here is restricted to staying in the
natural-like subspace; once you're in, the lever shrinks below noise.

## What this rules out
- Structured motif arrays as a way to teach grammar
- High-density motif training sequences (despite their high oracle scores)
- The hypothesis that the plateau was a GC-50-random-specific ceiling

## What to try next
Now that we know the plateau is robust and natural-like wins:
1. **Map composition curve** (016): Per-sequence GC variance test — does
   a bimodal library (half GC-30, half GC-70) hurt despite library mean=0.5?
2. **Multi-seed oracle pool** (017): Combine top selections from 5
   different seed runs of 007's recipe — does seed averaging reduce
   noise and beat single-seed?
3. **Cell-type-minimum oracle** (018): Select top by MIN(K562, HepG2,
   SKNSH) — sequences active in ALL three cells (pan-active).
   This biases for transferable activity.
4. **cCRE unfiltered** (019): Natural cCRE sequences with their
   intrinsic GC distribution — does natural composition variance help?
5. **Adversarial low-bound** (020): GC=0.30 single-direction —
   confirms downward boundary for symmetry analysis.
