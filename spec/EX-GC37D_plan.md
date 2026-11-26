# EX-GC37D — waves under realised random coincidences

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread ("close 37A/B. Go ahead with the next steps", 2026-11-26 01:31 UTC). Written before any run. Default cap 2 core-hours; expected under 0.3. This is our construction, not LaViolette's Model G.

## Question

37B's wave rule applied the *mean* coincidence weight N_xN_y on every face at every event. In the token picture, a transfer either happens or does not. Does a wave whose face transfers are realised stochastically still bend twice as much as lapse-only? And how much coherence does it lose?

## Method

- 37B B1b setup exactly: the 35K packet, mass source and readouts; h = .16; proper step .04.
- At each site event, each face transfer is **accepted** with probability p = p0·N_xN_y. An accepted transfer carries weight 1/p0, so the mean weight equals 37B's N_xN_y.
- p0 values: 1, 0.5, 0.25. The stability limit is a weight of 8 or less. Three seeds for each p0.
- **Controls:**
  - the same stochastic rule with lapse-only weights (accept with probability p0·N_f/N_x, weight 1/p0), with the same seeds;
  - flat N.
- **Readouts:**
  - group bend for each run;
  - the ratio of handshake to lapse-only bend, per seed and as a mean;
  - norm ratio;
  - coherent fraction, defined as the power in the forward carrier peak (|k − k₀| < 1) over total power;
  - the flat run's spurious angle.

## Predictions (scored)

- **P1:** the mean bend ratio is within 3% of 37B's deterministic 2.049 for every p0.
- **P2:** the flat-run spurious angle is under 0.05° for every p0 and seed, so noise does not mimic bending.
- **Descriptive:** coherent fraction against p0.

## Outputs

`src/exgc37d.py`, `results/EX-GC37D/`. The separate check: p0 = 1 with N ≡ 1 accepts every transfer, so it must reproduce the deterministic flat run bit for bit.
