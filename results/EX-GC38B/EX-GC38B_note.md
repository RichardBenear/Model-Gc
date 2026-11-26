# EX-GC38B — is the entanglement temperature the clock? (lattice Bisognano–Wichmann and Tolman)

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC38B_plan.md) (implementation spec A1 committed before any run). PI closure pending. Cloud container, Python 3.13, python-flint (arbitrary precision). About 0.4 core-h. This is our construction, not LaViolette's Model G.

**Bisognano–Wichmann holds on the lattice almost exactly. The Tolman response to a lapse has the right sign and shape but only about 60% of the predicted size.**
- **Flat lattice, half-chain cut.** The vacuum looks thermal with respect to the boost: the entanglement Hamiltonian is K ≈ 2π Σ x·(local energy density). The local inverse temperature is β(x)/(2πx) = 0.9998 at the cut, and still 0.992 twenty sites out.
- **With a static lapse bump N = 1 + 0.1f(x).** The local β rises where N is larger, i.e. the local temperature falls where clocks run faster, as Tolman requires. But the change is only 0.61–0.63 of the optical-length prediction. The RMS miss is 38%, against a 10% gate: **P2 FAIL.**
- The miss does not depend on the mass, and it does not shrink when the bump is made 2× or 3× smoother in lattice units (40%, 39%, 39%). Part of the response sits in non-local terms of K.

## Results (open chain of 200 sites, Dirichlet walls, cut in the middle; m = 0.05; 1000-bit arithmetic)

**P1, lattice Bisognano–Wichmann: PASS.** β_π(j) = K's π_j² coefficient (×2), against 2π(j − ½):

| j | 1 | 3 | 5 | 10 | 15 | 20 | 25 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| β/(2πx) | 0.99984 | 0.99983 | 0.99973 | 0.99852 | 0.99511 | 0.99161 | 0.98583 |

The gate was within 10% for 3 ≤ j ≤ 20; the range is 0.9916–0.9998. The slow fall-off is the walls and the mass (ξ = 20). The entanglement spectrum is evenly spaced (ε = 1.94, 5.83, 9.72, 13.61, …, spacing 3.89), as the corner-transfer-matrix result for this chain predicts.

**P2, Tolman with a lapse: FAIL.** f = exp(−(x − 10)²/32) on the region side, ε = 0.1. Prediction (A1): β_loc = 2πN(x)∫₀ˣdx′/N, so δβ = 2πε[x f − ∫₀ˣ f].

| x | 2.5 | 5.5 | 8.5 | 10.5 | 13.5 | 16.5 | 19.5 | 24.5 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| lattice δβ | 0.097 | 0.84 | 2.20 | 2.50 | 0.86 | −1.74 | −3.01 | −2.96 |
| predicted | 0.118 | 1.05 | 2.79 | 3.12 | 0.73 | −3.16 | −5.48 | −6.24 |
| naive β ∝ N | 0.27 | 1.84 | 4.98 | 6.55 | 5.78 | 2.77 | 0.73 | 0.02 |

- The RMS difference over 3 ≤ j ≤ 20 is **0.38** of the predicted change. The gate was 0.10.
- The naive "β ∝ N(x) at fixed x" form does much worse (0.87). The lattice does carry the non-local "optical distance" tail: it is negative beyond the bump, where N = 1 again. But it carries only about half of it (−3.0 against −6.2).

## Diagnostics (descriptive, run after P2 failed)

- **Not the mass.** The miss is the same for m = 0.05, 0.02 and 0.005 (0.375–0.377). It is also the same when the mass term is made uniform in optical length (m²/N instead of m²N), which removes the only term that breaks the conformal map.
- **Not the lattice scale.** With the bump scaled by s (centre 10s, width 4s, m = 0.1/s, L = 40s):

  | s | 1 | 2 | 3 |
  |---:|---:|---:|---:|
  | RMS miss | 0.405 | 0.391 | 0.394 |
  | best-fit amplitude | 0.613 | 0.629 | 0.626 |

  The amplitude settles near 0.63, so the miss is not a cutoff artefact that refinement removes.
- **A sharp step** N = 1 → 1.1 at x = 5 (L = 40) gives δβ = 2.2 just past the step, against 3.14 predicted. It then decays to 1.3 by x = 14.5, where the prediction stays constant. This holds for m = 0.005 and 0.1.
- **Part of the response is non-local.** The off-diagonal weight of K rises from 0.2% (flat) to 0.9% with the lapse. Summing K's π coefficients along a row brings the tail closer to the prediction (−4.7 against −5.5 at x = 19.5), but not near the bump. The π² and gradient readouts also differ from each other there. With a lapse, the lattice K is not of the strict form ∫β(x)ℋ(x).

## Checks

- **G1: PASS.** The entropy from K's spectrum (arbitrary precision) is 0.499664340328. The double-precision symplectic entropy of the same block is 0.499664340328, a difference of 1e-13.
- **G2, approach to the continuum: PASS.** At fixed x·m = 0.5, β/(2πx) is 0.99866, 0.99831 and 0.99777 for m = 0.1, 0.05 and 0.025. These are all within 0.25% of 1. The small drift comes from the walls (L = 50, 100 and 200 sites, i.e. 5ξ).
- **Precision:** 1000 and 1400 bits agree to 1.3e-6 in β. 500 bits, as first specified in A1, was too few: test runs before scoring showed the profile breaking down beyond x ≈ 10. So 1000 bits was used. This is a deviation from A1.
- **Consistency:** a uniform lapse on the whole chain leaves K unchanged, to all printed digits.

## What it means (inference)

- **The vacuum's entanglement does define a local temperature and a local "time": the boost.** On our lattice this holds almost exactly, right down to one site from the cut. So the clock that entanglement defines near any surface is the Rindler (Unruh) clock. This is the ingredient Jacobson's argument needs: the local first law δS = δ⟨K⟩, with K the boost energy.
- **The link to the lapse is only partial.** A lapse changes the entanglement temperature in the Tolman direction, and the change remembers the optical distance to the cut. But on the lattice it falls short by about 40%. Either the continuum formula needs a correction I have not identified, or the lattice modular Hamiltonian spreads the response into non-local terms that a local readout does not see.
- **For the gravity question:** 38A and 38B together say that the vacuum's entanglement gives the infrared part of gravity's stiffness and a boost clock at every surface. Neither shows how the lapse (37B's firing rate) is *set* by matter. That is the Jacobson step: use δS = δ⟨K⟩ for every small cut to derive the field equation for N. It is the natural next design note, not yet a run.

## Files

[Code](../../src/exgc38b.py) · [scaling diagnostic](../../src/exgc38b_scale.py) · [row-sum diagnostic](../../src/exgc38b_rowsum.py) · `main.json` (P1, P2, G1, precision), `g2.json`, `p2diag.json`, `scale_s1..3.json`, `rowsum.json`
