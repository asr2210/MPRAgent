# Exp 029: mixture of 10k each from seeds {1, 7, 42, 100, 2024}
- eval_01=**0.5194** K562=0.9945 HepG2=0.5640 SKNSH=-0.0004

Mixture roughly averages individual seed scores. As expected — each seed is i.i.d.
random uniform, so mixing them is statistically identical to a fresh single seed.
Useless variance reduction strategy in this regime.

(Individual seeds: 0.521, 0.518, 0.522, 0.520, 0.520 → mean 0.520, mixture 0.519)
