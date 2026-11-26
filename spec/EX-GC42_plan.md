# EX-GC42 — a pair network coupled to an energy-sourced lapse: does a hub source a universal 1/r?

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread: "Go ahead build: couple a pair network to a lapse field sourced by its energy density, and test whether a hub sources a universal 1/r lapse" (2026-11-26 17:27 UTC). Written before any run. Cap 4 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why

- **40C** found that a hub in a pair network is a local defect: every response decays exponentially. A long-range field needs a gapless carrier coupled to energy.
- **37C and 39C** require that carrier, the lapse N, to be relativistic but constrained by matter, not a free radiating field.

Here the two are joined. The lapse is fixed by a constraint equation sourced by the network's energy density, and the network feels the lapse through 37B's event weights.

**Honest framing.** If the source is a conserved energy density and the constraint is Poisson's equation, then 1/r with coefficient GΔE₀ follows at first order in G by construction. The test is therefore in what can spoil it:
1. **Does the vacuum gravitate?** The pair vacuum's own energy responds to N. If that response feeds back into the source, it acts like a mass term (or an anti-mass term) for N, and the 1/r law is destroyed beyond some length.
2. **Is the coefficient universal?** The vacuum's non-local response renormalises G. The question is whether that renormalisation is the same for every hub.
3. **Does a hub obey the equivalence principle under these weights?** That is, is its proper energy the same at every uniform N? 39B found this depends on how matter is built.

## Model

**Network.** 40A's pair network: periodic cubic lattice, six fermion modes per site, half filling. The background on-site coupling is λ₀ = 0.1 (main) and 0.2 (repeat).

**Coupling to the lapse (37B event weights).**
- Each link's hopping is weighted N_v N_w: a two-site coincidence.
- Each on-site term, including the hub's extra Δλ, is weighted N_v: a one-site event.
- So H[N] = −Σ_links N_v N_w (a† a + h.c.) + Σ_v N_v (λ₀ + Δλ_v) H₁,v.

**Energy density.** ρ_v[N] = Re Σ_{i∈v} (h[N] P)_ii, the on-site term plus half of each adjacent link, as in 40A. Σ_v ρ_v = E exactly.

**Constraint.** Lattice Poisson, Δn = 4πG S, with n = N − 1. This is a constraint, not a wave equation (39C). The source S comes in two variants, both fixed now:
- **Variant B (raw):** S_v = ρ_v[N] − ρ_vac[N ≡ 1]. All energy gravitates, including the vacuum's response to N.
- **Variant A (local vacuum subtracted):** S_v = ρ_v[N] − ρ_vac(N_v). Here ρ_vac(N) is the vacuum energy density at uniform lapse N. This subtracts the vacuum's local response, which is the lattice version of tuning the cosmological constant. It is put in by hand, and the note will say so.

**Hubs.** One site, or the 7-site cluster, with Δλ ∈ {0.1, 0.2, 0.4}. Coupling strengths 4πG ∈ {0.01, 0.05}.

## Method

- **42A (uniform N).** Bloch H at N = 0.9, 1.0 and 1.1. Record ρ_vac(N), and the hub's ΔE₀(N) and bound-state energies from 40C's Woodbury method. Proper quantities are those divided by N, the local clock rate.
- **42B (vacuum response kernel K(q)).**
  - Linear response of S to a lapse perturbation n_q: δS_q = K(q) n_q.
  - It is computed exactly with a supercell. Set n(x) = ε cos(qx) along x, with Bloch momenta in y and z, at L = 32 and q = 2πm/32, m = 1..8. Use ε = ±1e-4 and central differences.
  - K(0) comes from the uniform-N derivative.
  - One direction check: q along (1,1,0), using an L = 16 supercell in the x+y direction, for m = 1, 2.
- **42C (hub far field, linear order in the coupled system).**
  - n_q = −4πG S_hub,q / (q̃² + 4πG K̃(q)), with q̃² the lattice Laplacian eigenvalue, and K̃ = K (variant B) or K − K(0) (variant A).
  - S_hub(x) is the hub's flat-space energy-density change from 40C's Woodbury method, taken from the site blocks of δP (one site and 7 sites).
  - n(x) on a 64³ grid by FFT, with K̃(q) interpolated in |q| from 42B.
  - Fits outside the core (r = 4–20):
    - Variant A: n = −G_eff M/r against the periodic lattice Green's function.
    - Variant B: a damped or oscillating form.
- **Check C2 (does linear order hold?).** Dense self-consistent solve at L = 8, alternating "diagonalise H[N], then solve Poisson" with mixing 0.5, for one hub at 4πG = 0.01, in both variants. Compare n at r = 2, 3, 4 with the linear-order prediction computed on the same L = 8 grid. For that, the exact K kernel on L = 8 comes from a dense one-site lapse perturbation.

## Predictions

- **P1 (equivalence principle under 37B weights).** Expectation: **fails.** The hub's proper energy ΔE₀(N)/N changes by more than 1% between N = 0.9 and 1.1.
  - Reason: hopping scales as N² and the on-site term as N, so the ratio of hopping to λ changes with N, and the proper physics changes. A lattice-scale hub is not "built like light" in 39B's sense.
  - Recorded: the exponent of ΔE₀/N in N.
- **P2 (the raw vacuum gravitates, variant B).** K(0) = dρ_vac/dN ≠ 0, with K(0) < 0 at λ₀ = 0.1 and 0.2.
  - Reason: by Hellmann–Feynman, dρ_vac/dN = 2e_link + e_site, and both are negative (bonding energies).
  - Then the variant-B far field is oscillatory, cos(k_J r)/r with k_J² = 4πG|K(0)|, not 1/r. Scored: the fitted k_J within 10% of √(4πG|K(0)|) at both G.
  - Meaning: raw vacuum energy gravitating destroys Newton's law. This is the cosmological-constant problem in miniature.
- **P3 (variant A gives universal 1/r).** For every hub (1 and 7 sites, Δλ = 0.1, 0.2, 0.4), both λ₀ and both G:
  - n is 1/r outside the core: the periodic-Green's-function fit residual is below 1% of n at r = 4–20;
  - and G_eff M_fit / (G ΔE₀) is the same for all hubs, within 1%, at fixed λ₀ and G. This is universality: one G_eff, with every hub's mass equal to its energy.
  - G_eff/G = 1/(1 + 4πG c₂), with c₂ from K(q) − K(0) ≈ c₂ q̃². The sign of c₂ is recorded, not predicted.
- **P4 (descriptive).** The size of the G renormalisation, and whether K(q) − K(0) is analytic in q² (gapped network, expected) or has a non-analytic part.

## Checks

- **C1.** Σ_v S_hub,v = ΔE₀ (energy partition), to 1e-8 at L = 16.
- **C2.** Dense self-consistent vs linear order, at L = 8 and 4πG = 0.01. The values of n at r = 2–4 agree within 5% (the O(G²) corrections should be about 1%).
- **C3.** The supercell K(q) at q → 0 tends to the uniform-N derivative K(0), within 2% at m = 1 after the q² trend is extrapolated out. K along (1,1,0) agrees with K along x at equal |q̃| within 5% at small q (isotropy).

## Outputs

`src/exgc42.py`, `results/EX-GC42/`.
