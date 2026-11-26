# EX-GC39B — does the lapse shape matter? (equivalence-principle scaling; self-held lapse wells)

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC39B_plan.md) · [design note](../../spec/EX-GC39_design.md). PI closure pending. Cloud container, Python 3.13. About 0.2 core-h. This is our construction, not LaViolette's Model G.

**The lapse can shape matter consistently, and can even hold matter together on its own. But only matter built like light does it right, and the lapse does not pick out particle masses.**
- **B1: what kind of matter passes.** Under 37B's event weights, particles made of waves keep their proper size and clock in any lapse and fall at GR's rate. A φ⁴ kink falls at 0.997 of the GR rate, and two different kinks agree to 0.03%. Reaction–diffusion particles (Model Gc) do not: their proper size scales as N^{-1/2}, so local physics changes from place to place. The Skyrmion passes only if its four-derivative term is carried by events of weight N⁶.
- **B2: can the lapse hold a particle together?** Yes, it can bind a field into a stable lump with no other binding force. But the lumps form a continuous family (the Schrödinger–Newton/boson-star family), and it ends at a maximum mass GμM ≈ 0.55–0.6, where the lapse collapses. Apart from that endpoint, no mass is selected.

## B1a — Model Gc particle (3D radial, k = 0.1, G feeds X)

Reactions are one-site events (rate N); diffusion hops are coincidences (rate N²).

| N | W rms radius, coordinate | proper (÷N) | slow eigenvalue, proper |
|---:|---:|---:|---:|
| 0.8 | 1.2822 | 1.6028 | −0.00221 |
| 0.9 | 1.3601 | 1.5112 | −0.00249 |
| 1.0 | 1.4337 | 1.4337 | −0.00276 |
| 1.1 | 1.5037 | 1.3670 | −0.00304 |

- **P1 (the expected failure of the equivalence principle): confirmed.** The proper radius scales as N^{−0.4996}. A Gc particle deeper in a lapse well is larger in its own units, so it would not fall universally.
- **P2: FAIL.** The slowest relaxation rate spreads by 32% instead of staying within 1%. The leading eigenvalue is a slow mode of the carrier diffusing out to the fixed outer boundary (R = 60), with rate ∝ D_eff/R² ∝ N, not an internal clock of the particle. That is the same defect seen another way: diffusive dynamics do not scale with the lapse.
- **Check:** at N = 1 the solver reproduces the ex008 particle exactly (core W 3.15071).

## B1b — φ⁴ kink, 1+1D lattice (h = 0.1)

| Rule | Proper-width exponent | Proper shape-mode frequency exponent | Fall, a/(−g), kink 1 | kink 2 |
|---|---:|---:|---:|---:|
| Handshake (gradient weight N_iN_j, potential N²) | −0.002 / +0.006 | 0.0003 / 0.0006 | 0.9976 | 0.9973 |
| Lapse-only (potential weight N) | +0.50 / +0.49 | −0.50 / −0.50 | 0.4997 | 0.4999 |

Kink 1 is (λ, v) = (1, 1); kink 2 is (8, 0.5), with a different mass and width. The "/" pairs give the two kinks.

- **P3: PASS.** With the covariant rule, proper width and frequency are unchanged to 0.6% or better. Both kinks fall at 0.997 of the geodesic rate and agree with each other to 0.03%.
- **P4: half FAIL.** The lapse-only rule does break proper size, as predicted (∝ N^{1/2}). But it falls *universally*, at exactly half the rate, rather than non-universally as I predicted. It behaves as if the lapse were √N. The registered "two kinks differ by more than 5%" does not hold: they agree to 0.04%.
- **Checks:** in flat space the kink stays put (drift 3e-13) and energy is conserved to 1e-9.

## B1c — 3D Skyrmion hedgehog (radial; pion mass 0.5, Skyrme units)

| Quartic-term weight | proper-size exponent | proper-energy exponent |
|---|---:|---:|
| N⁶ (covariant, p = 6) | 5e-5 | 5e-5 |
| N⁴ (two coincidences, p = 4) | −0.82 | −1.08 |

- **P5: PASS.** With p = 6 the Skyrmion is identical in proper units at N = 0.9, 1.0 and 1.1 (radius 0.91402, energy 154.21). With p = 4 it shrinks and gets lighter in proper units in deeper wells.
- Grid check: E = 154.2102 at n = 600 and 154.2121 at n = 1200 (1e-5).

## B2 — a lapse well holding its own matter

Complex scalar of mass μ in the metric −N²dt² + N⁻²dx², with N sourced by the scalar's energy (∇²N = 4πGρ). Radial lattice, h = 1, R = 4000. C = Gμ²M·R90.

| GμM | 0.08 | 0.16 | 0.3 | 0.5 | 0.55 | ≥ 0.6 |
|---|---:|---:|---:|---:|---:|---|
| C (μ = 0.1) | 6.80 | 6.68 | 6.26 | 4.89 | 4.14 | no solution: the lapse collapses |
| N(0) | 0.998 | 0.992 | 0.970 | 0.887 | 0.836 | → 0 |
| ω/μ | 0.999 | 0.996 | 0.985 | 0.944 | 0.919 | |

- **P6: PASS.** In the Newtonian regime (1 − N(0) < 0.05, R90 ≤ 1000), C = 6.82–6.26 for all μ (0.05, 0.1, 0.2) and for G = 21. That is constant within 8%, against a 10% gate. A separate continuum Schrödinger–Newton solve gives C = 6.833, and the weakest lattice lumps give 6.81–6.82. The lattice cutoff plays no role while R90 ≫ h.
- **P7 (descriptive): the family has a maximum mass.** Self-consistent lumps exist up to GμM ≈ 0.55 and fail between 0.55 and 0.6 (finer steps, slower mixing). Beyond that the central lapse is driven to 0. This is the analogue of the boson-star maximum (0.633 in GR). The lumps are node-free ground states; the first radial excitation also exists (C = 32.5 at GμM = 0.3).
- The masses are not selected. For any μ and G, any M below the maximum gives a lump, with size R90 ≈ 6.8/(Gμ²M). With the lattice G (≈ 21 from 37A), a lump of mass μ (one quantum) would be about 6.8/(21μ³) lattice units across, i.e. enormous for light fields.
- **Checks:** G = 0 gives no bound state (ω/μ = 1.0001). G = 21 at fixed GμM reproduces G = 1 exactly (the expected scaling).

## What it means (inference)

- **On your idea that the lapse shapes how matter emerges:**
  - It is consistent only if matter is made of the same kind of events as light. Its internal "springs" (gradient terms) must be two-site coincidences, its on-site terms one-site events, and its higher-derivative terms must carry the matching higher power of N.
  - Then the lapse shapes matter invisibly to local observers (same proper size and clock everywhere), and matter falls along geodesics automatically. That is the equivalence principle emerging from event counting, not imposed.
  - Model Gc's reaction–diffusion particles fail this. A future particle on this substrate should be wave-built (Skyrmion-like), as GC26 already suggested.
- **The lapse can hold matter together.** A self-sourced lapse well traps its own field and makes a stable particle-like lump. But the lump's mass is free, so the lapse does not explain particle masses. The only mass it fixes is a maximum, M_max ≈ 0.55/(Gμ).

## Files

[B1 code](../../src/exgc39b1.py) · [B2 code](../../src/exgc39b2.py) · `b1a.json`, `b1b.json`, `b1c.json`, `b2.json`, `b2_edge.json`
