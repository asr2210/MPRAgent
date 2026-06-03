# 008_reporter_construct

Mimics typical synthetic MPRA library: 140bp enhancer with 4-8 TF binding
motifs in random background + 60bp fixed minP cassette (TATA-Inr-spacer).

Hypothesis: providing the stereotyped construct architecture that most MPRA
training data has might let v14's model latch onto a useful feature.

Result: eval_01 = -0.0039. Noise. Identical to v14's random_uniform baseline.

**Conclusion**: structural priors about reporter architecture do not help.
