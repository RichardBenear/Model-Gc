# EX-GC44 — separating the clock from space: lapse homogeneity, a spatial link factor, and whether a bound composite is stress-free

2026-11-26 | EXPLORATORY | Claude | Plan only. PI request in thread: "Give me you plan based on that review" (2026-11-26 23:54 UTC). **Not run. Running needs a separate PI go-ahead.** It follows the [Codex review of GC43](../review/EX-GC43_codex_review.md). Cap 2 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why

43 found Q = ΔE + ΔE_link for the 37B weights: links carry N², sites carry N. Codex's review makes three points:
1. Turn the equal-energy comparison into an end-to-end far-field benchmark.
2. Use a uniform-lapse scaling control: if H(N) = N·H(1), then Q = ΔE identically. Then ask whether a *separate spatial rule* can keep the 41 optics.
3. Do not invoke a von Laue cancellation without a stress observable defined from strain, independently of the lapse.

**Reading (inference, written before any run).**
- The three points fit one structure, the one GR itself has. In static weak-field GR there are two potentials:
  - Ψ (space), sourced by energy: ∇²Ψ = 4πGρ;
  - Φ (time), sourced by ρ + 3p: ∇²Φ = 4πG(ρ + 3p).
- Light feels both, through Φ + Ψ. Massive bound bodies have ∫p = 0 (von Laue), so the two potentials agree.
- 37B effectively set Ψ ≡ Φ by putting N² on links. 43's charge, ΔE + ΔE_link, then looks like ρ + 3p, with ΔE_link in the role of 3p.
- If this is right, 43 is not a contradiction of GR. It says this network's matter is radiation-like: its energy, even at rest, sits largely in hopping.
- The test is whether the split survives when it is built explicitly, and whether any bound composite has zero lattice stress.

## Model

- **Network.** The 40A pair network, as in 42 and 43, at λ₀ = 0.1 (main) and 0.2 (repeat).
- **Two fields per lattice.** A lapse N_v on sites (the clock rate) and a spatial link factor a_vw on links (the proper length of a link in coordinate units).
- **Matter weights (44B onward):**
  - H[N, a] = Σ_v N_v H_site,v + Σ_links (N̄_vw / a_vw) H_link,vw, with N̄_vw = (N_v + N_w)/2.
  - This is homogeneous of degree 1 in N, so under uniform N the whole Hamiltonian scales as N.
  - The signal (hop) rate per coordinate link is N̄/a, so light's coordinate speed is N/a.
- **Charges.**
  - Lapse: Q_N = Σ_v ∂H/∂N_v, which equals ΔE identically at N = a = 1.
  - Spatial factor: Q_a = Σ_links −a ∂H/∂a = ΔE_link. This is the "stress" defined from strain: Codex's independent observable, with hop ∝ 1/a as the stated definition.
- **Field equations.** These are linearised static GR, with Φ = n = N − 1 and Ψ = −ln a:
  - ∇²Ψ = 4πG ρ_N
  - ∇²Φ = 4πG (ρ_N + ρ_a)
  - The densities ρ_N and ρ_a are the site densities of Q_N and Q_a. The local vacuum counterterm (43's prescription L) is applied to both.
  - These equations are put in; they are not derived from the network. The test is the matter side: what the two charges are, and whether they can be equal.
- **Mapping check.** With a ≡ 1/N̄ imposed (Ψ ≡ Φ), the link weight becomes N̄², and 43's charge ΔE + ΔE_link is recovered exactly (check C2).

## Parts

### 44A — end-to-end equal-energy benchmark (Codex item 1, 43's model unchanged)

- **Sources.** Neutral particle–hole pairs at fixed particle number:
  - Pair P2: top (dispersive) band packet, σ = 2, × top valence band packet, σ = 2.
  - Pair P1: flat conduction band × top valence band packet. The valence σ is tuned by bisection so that ΔE(P1) = ΔE(P2) to 1e-10. The flat band's energy is fixed, so only the hole can be tuned.
  - The existing charged-vs-neutral coincidence at λ₀ = 0.2: top-band packet vs particle–hole, both at ΔE = 0.8791069.
- **Preparation.** Static t = 0 states. Q is exactly conserved in time for these states: each band packet has no interband terms. The 200-time check of 43's C3 is repeated for one pair.
- **Profiles.** The spatial conjugate-charge profile q_v(t = 0), from real-space ⟨∂H/∂N_v⟩ at L = 32, embedded in a 128³ periodic grid. The lapse is solved under 43's prescription L, at 4πG = 0.01.
- **Readout.** The far-field amplitude ratio between sources, from the periodic-Green's-function fit at r = 6–20.

### 44B — uniform-lapse scaling control (Codex item 2)

- Matter weights N̄ on links, N on sites, with a ≡ 1.
- Recompute Q/ΔE for all 43 packets, pairs and hubs.
- Recompute ΔE₀(N)/N for the hubs at N = 0.9 / 1.0 / 1.1.
- Optics: the group velocity of a top-band packet at uniform N = 0.9 / 1.0 / 1.1, fitted to v ∝ N^p.

### 44C — the two-field split

For every 43 state and hub, at N = a = 1:
- Q_N and Q_a, from finite differences of H in uniform N and uniform a, separately.
- The far-field coefficients of Ψ (∝ ΣQ_N) and Φ (∝ ΣQ_N + ΣQ_a).
- The light-bending coefficient Φ + Ψ, relative to a "dust" source of the same energy (Q_a = 0).

### 44D — a bound composite: is its excess stress zero? (Codex item 3)

- **The composite.** An exciton of the unchanged network: one particle in the top (dispersive) band and one hole in the top valence band.
  - The interaction is an on-site density–density attraction, V = −U Σ_v n_v^particle n_v^hole.
  - It is solved exactly in the band-projected two-body problem at total momentum K = 0, at L = 16 (4096 relative momenta, dense diagonalisation).
- **Mobility.** The exciton's dispersion at K = 2π/16·(1,0,0) and (2,0,0) is recorded, to show it is mobile: Codex's "mobile bound state".
- **Strain response of the interaction.** U ∝ a^{−w}, scanned over w ∈ {0, 0.5, 1, 2}.
- **Couplings.** U is chosen by bisection so that the binding energy E_b is 5%, 15% and 30% of the pair gap.
- **Measured for each:**
  - ΔE (above the vacuum);
  - ΔE_link, by Hellmann–Feynman (a finite difference in the hopping scale);
  - ⟨V⟩;
  - the excess stress δQ_a(w) = Q_a(bound) − Q_a(unbound pair at the same band edges) = δ⟨H_link⟩ + w δ⟨V⟩;
  - the "von Laue exponent" w* at which δQ_a = 0.

## Predictions (scored)

- **P1 (44A).**
  - Far-field amplitude ratios equal the charge ratios within 1e-3. The far field is pure monopole in Q at linear order, so this checks the end-to-end chain.
  - The equal-energy neutral pairs P1 and P2 differ in amplitude by more than 5%.
  - The charged-vs-neutral case reproduces the 56% difference within 1%.
- **P2 (44B).**
  - Q/ΔE = 1 for every state and hub, within 1e-8.
  - ΔE₀/N changes by less than 1e-8 between N = 0.9 and 1.1: the equivalence-principle redshift test passes.
  - The optics exponent p = 1.000 ± 0.005, so light has only the one-site rule's Newtonian (factor 1) bending.
- **P3 (44C).**
  - Q_N = ΔE (within 1e-8) and Q_a = ΔE_link (within 1e-6) for all states.
  - The Φ + Ψ bending coefficient relative to dust equals 1 + ΔE_link/(2ΔE), between 1.4 and 1.7 for single packets.
  - Light from this matter bends GR-correctly only if ΔE_link = 0, which no 43 state satisfies.
- **P4 (44D, the main question).**
  - A bound exciton exists for U above a threshold, and is mobile: its dispersion is nonzero at K ≠ 0.
  - There is no lattice von Laue identity: w* changes by more than 10% between E_b = 5% and 30% of the gap.
  - If w* is constant within 2%, that is a surprise, and the note says so: it would mean a strain scaling for the interaction exists that makes binding stress-free.
- **P5 (descriptive).** The constituents' own Q_a at the band edges: what fraction of each constituent's rest energy is link energy. This is what makes the matter radiation-like, before any binding.

## Checks

- **C1.** The real-space profile sums ΣQ_v equal the band formula of 43 to within 1e-10, for each 44A source.
- **C2 (mapping).** With a ≡ 1/N̄ imposed, the uniform-N derivative of the 44 Hamiltonian equals 43's Q for all packets and hubs, within 1e-8.
- **C3 (two-body solver).** At U = 0 the spectrum equals the sums e_particle − e_hole. For one bound state at L = 12, the band-projected energy is compared with a direct real-space Lanczos solve of the same projected problem, within 1e-8.
- **C4 (Hellmann–Feynman).** ΔE_link from the hopping-scale derivative of E_bound equals ⟨ψ|∂H/∂t|ψ⟩ from the eigenvector, within 1e-8.
- **C0 (field equations).** The two linear equations are solved for a unit point source of each charge. The far-field coefficients of Ψ and Φ are checked within 1e-3 of the analytic values.

## What would change the picture

- **If P4 holds** (no lattice von Laue), a fixed-lattice network with lapse-weighted hopping cannot make massive bodies whose active mass equals their energy. The space–time split is right, but the matter is radiation-like and binding does not cure it. The next substrate would need matter whose rest energy is on-site (Q_a = 0 at rest), with binding governed by a dilation-like symmetry.
- **If some w* is universal,** that fixes how interactions must respond to spatial strain. That would be a microscopic derivation of the stress cancellation, and the first concrete mechanism-level statement this programme has produced about why massive bodies gravitate by their energy.

## Outputs

- `src/exgc44.py`, `results/EX-GC44/`.
- A separate recomputation, `src/exgc44_recompute.py`, recomputes Q_N, Q_a and w* from eigenvalue finite differences only, with no eigenvectors.
