# EX-GC39B — does the lapse shape matter? (B1 equivalence-principle scaling; B2 self-held lapse wells)

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread (2026-11-26 04:53 UTC). Written before any run. Cap 4 core-h. Design: [EX-GC39_design.md](EX-GC39_design.md). This is our construction, not LaViolette's Model G.

## B1 — do emergent particles keep their proper size and clock in a lapse?

37B's rule: one-site events run at rate N, and two-site coincidences (hops, propagation) at rate N². The weak-field metric it gives is −N²dt² + N⁻²dx², so proper length is dl = dx/N. The equivalence principle needs coordinate size ∝ N, proper size and proper frequency independent of N, and universal fall.

**B1a, Gc1 reaction–diffusion particle** (`ex008_longrange.Model`, 3D radial, k = 0.1, G feeds X).
- In uniform N, reactions run at rate N and diffusion (hops) at N², so in the steady state every diffusion constant is multiplied by N.
- Solve for N = 0.8, 0.9, 1.0, 1.1. Measure the W rms radius (coordinate), the proper radius = coordinate/N, and the leading l = 0 relaxation eigenvalue (coordinate rate; proper rate = rate/N).
- **P1:** proper radius ∝ N^{-1/2} (fitted exponent −0.5 ± 0.05). The equivalence principle FAILS.
- **P2:** proper relaxation rate independent of N within 1%.

**B1b, a wave-built (kink) particle, 1+1D lattice, with falling.**
- φ⁴ kink, V = (λ/4)(φ² − v²)². The covariant event rule (handshake): φ_tt(i) = N_i²[Σ_j N_iN_j(φ_j − φ_i)/h² − V′(φ_i)]. This is the discretisation of −N²dt² + N⁻²dx², with face weight N_iN_j as in 37B.
- Control rule ("lapse-only"): potential term at rate N_i (weight N_i instead of N_i²).
- h = 0.1, leapfrog dt = 0.02.
- Uniform N = 0.8, 0.9, 1.0, 1.1: measure the static width (coordinate) and the small-oscillation (breathing/shape-mode) frequency.
- Fall: N = 1 + g·x with g = 1e-3, two kinks with (λ, v) = (1, 1) and (4, 0.5), i.e. different masses and widths, released at rest at x = 0. Measure the coordinate acceleration from a quadratic fit of the centre (zero crossing) over t ≤ 40. The geodesic prediction is a = −N³∂_xN ≈ −g.
- **P3 (handshake):** proper width and proper shape-mode frequency independent of N within 1%; both kinks fall at a = −g within 2%, and agree with each other within 1%.
- **P4 (lapse-only control):** width ∝ N^{+1/2} and non-universal fall (the two kinks differ by more than 5%).

**B1c, the 3D Skyrmion: which event weight the quartic term needs.**
- Radial hedgehog energy E = ∫[w2·(sigma-model term) + w4·(Skyrme term) + w0·(pion-mass term)] with weights from the event rule: w2 = 1, w0 = N⁻², and w4 = N^{p−4} for a quartic term carried by events of total weight N^p. Covariance needs p = 6 (w4 = N²).
- Minimise the radial profile for N = 0.9, 1.0, 1.1, for p = 4 and p = 6, with pion mass 0.5 (Skyrme units).
- **P5:** p = 6 gives proper size and proper energy·N⁻¹ constant within 1%. p = 4 does not; its size deviation is reported.

## B2 — can a lapse well hold its own matter?

- **Model.** A massive scalar χ, mass μ, in the static isotropic metric −N²dt² + N⁻²dx² (covariant action √−g = N⁻², weights N⁻⁴χ̇², (∇χ)², N⁻²μ²χ²).
- **Mode equation:** ∇²u + (ω²N⁻⁴ − μ²N⁻²)u = 0.
- **Lapse:** sourced by the mode's energy, ∇²N = 4πG ρ, with ρ = ½(ω²N⁻⁴u² + |∇u|² + μ²N⁻²u²) and N → 1 at large r (the 37B weak-field relation, extended as N itself).
- **Lattice:** radial, h = 1 (the cutoff), out to R = 4000, Dirichlet.
- **Self-consistency:** fix the total mass-energy M = ∫ρ d³x, iterate to convergence. Scan M at G = 1 and μ = 0.05, 0.1, 0.2, and also at G = 21 (37A's G_ind) with μ = 0.1.
- **Readouts:** existence and node-free ground state; ω; R90 (radius holding 90% of ρ); the central lapse N(0); and C = G μ² M R90 (in the Newtonian Schrödinger–Newton limit, C is a constant, about 9.9 for R99; R90 is computed from the converged profile).
- **P6:** lumps exist for all scanned M with N(0) > 0.5. C is constant within 10% while R90 > 20h and 1 − N(0) < 0.05. This is the Newtonian, continuous family: no mass is selected.
- **P7 (descriptive):** where the family ends. Either a maximum M (relativistic or cutoff), or R90 reaching a few h. Report M_max·Gμ and whether anything discrete appears besides the radial excitations (nodes).

## Checks
- **B1a:** N = 1 reproduces the ex008 particle (core W and consumption, to 1e-8).
- **B1b:** with N ≡ 1, the kink stays at rest (|x_c| < 1e-6 over t = 40) and energy is conserved to 1e-6.
- **B1c:** the p-independent N = 1 hedgehog energy matches the known massless value 1.232 × 12π² (Skyrme units, 6π²·… convention stated in code) within 1%, or the code's own fine-grid limit.
- **B2:** with G = 0, the mode equation has no bound state. In the weak limit, a separate shooting solve of the Schrödinger–Newton equations (continuum, scaled) gives the same C within 2%.

## Outputs
`src/exgc39b1.py`, `src/exgc39b2.py`, `results/EX-GC39B/`.
