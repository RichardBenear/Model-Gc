# EX-GC37D — waves under realised random coincidences

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC37D_plan.md). PI CLOSED 2026-11-26. Cloud container, Python 3.13, numba. About 0.3 core-h. This is our construction, not LaViolette's Model G.

**Fails.** If each face transfer of the wave is realised as a random yes/no event, the randomness pumps energy into the wave.
- **p0 = 0.5 and 0.25:** the field grows by 10¹⁷ and 10⁵⁰. That is noise-driven, parametric growth.
- **p0 = 1** (99% of handshake transfers accepted): the packet survives but grows by 40–55%, keeps only 39–52% of its power in the forward carrier peak, and bends 2.64–2.77° against 3.31° for the deterministic rule.

The 37B factor of two holds for the *mean* coincidence rate and for independent tokens. It does **not** hold for a coherent wave driven event by event.

## Results (37B B1b setup: h = .16, proper step .04; three seeds each)

| p0 | Rule | Group bend (3 seeds) | Norm ratio | Coherent fraction | Acceptance |
|---:|---|---|---:|---:|---:|
| 1 | flat | 0.000°, 0.000°, 0.000° | 1.0016 | 1.00 | 1.000 |
| 1 | lapse only | 1.618°, 1.617°, 1.622° | 0.975 | 0.995 | 1.000 |
| 1 | handshake | 2.766°, 2.736°, 2.640° | 1.39–1.55 | 0.39–0.52 | 0.992 |
| 0.5 | all three rules | blow-up | ~10¹⁷ | <0.01 | 0.50 |
| 0.25 | all three rules | blow-up | ~10⁵⁰ | <0.02 | 0.25 |

- **P1 (mean bend ratio within 3% of 2.049): FAIL.** At p0 = 1 the ratio is 1.68 (2.714/1.619). At smaller p0 there is no packet left to measure.
- **P2 (flat spurious angle < 0.05°): FAIL** at p0 = 0.5 and 0.25, where the run blew up. It holds at p0 = 1, but trivially: with N ≡ 1 every transfer is accepted.
- **Separate check: PASS.** The p0 = 1 flat run reproduces the deterministic 37B flat run exactly (norm 1.00157949300422 and angle −7.161e-7° in both).
- **Lapse-only at p0 = 1 is nearly deterministic.** N_f/N_x ≥ 1 on about half the faces, so the acceptance probability is capped at 1 there. Those faces carry weight 1 instead of N_f/N_x, a bias of order 1%.

## Why (inference)

- **Each accepted or refused transfer is a random kick to the wave's stiffness at that moment.** A leapfrog wave with randomly time-varying stiffness is not energy conserving. Its energy grows exponentially, faster for larger weight fluctuations (1/p0). This is noise heating of the numerical scheme, and also of any medium built this way.
- **The counting mechanism must therefore act on many events before it acts on a wave.**
  - The coherent field has to be a coarse-grained average over many coincidences, like the token ensemble in 37B B1a, where each token is independent.
  - It cannot be one amplitude updated by each individual coincidence.
- **This is a constraint on any microscopic version of 37B:** light would be a collective mode of the coincidence events, not a field kicked by them.

## Files

[Code](../../src/exgc37d.py) · results `flat_*.json`, `lapse_*.json`, `handshake_*.json` (per-run hashes included)
