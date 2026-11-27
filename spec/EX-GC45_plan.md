# EX-GC45 — heavy on-site matter: can a body's spatial charge be small, and does its binding energy gravitate as energy?

2026-11-27 | EXPLORATORY | Claude | Plan only. PI request in thread: "Read and determine a next test" (2026-11-27 01:31 UTC). PI run go-ahead 2026-11-27 01:37 UTC ("ok, go ahead with 45"), given after this plan was pushed and before any run. Follows the [Codex review of GC44](../review/EX-GC44_codex_review.md), items 1 and 2. Item 3, a two-field optics test, is deferred to a later study. Cap 2 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why

44 split the coupling into two charges.
- **Clock charge Q_N.** Under lapse-homogeneous weights it equals ΔE for every state.
- **Spatial (strain) charge Q_a.** This sources the time potential Φ = Q_N + Q_a. In the pair network, Q_a ≈ ΔE even at rest, because all matter energy is hopping.

Codex's reanalysis corrects 44's reading. Cancelling the *whole* exciton's Q_a needs positive, strongly state-dependent strain exponents. Negative exponents cancel only the binding increment.

So the useful questions are two.
1. Can matter whose rest energy is on-site have small Q_a/ΔE?
2. For a bound body of such matter, does the binding energy carry its own Q_a?

GR's answer to 2 is no: a bound body's total stress is zero. A residual δQ_a ≈ η E_b would be a lattice Nordtvedt effect: bodies with different binding fractions fall differently, at order η·E_b/ΔE.

## Model

**Matter.** The 40A network plus one heavy mode per site, in two species (s = 1, 2, like spin). Each species is an independent copy:
- 6 link modes with hopping ∝ N̄/a, as in 44;
- on-site λ₀(J − I) ∝ N;
- a heavy mode h_v at energy m₀ ∝ N;
- an on-site hybridisation g between h_v and each of the site's 6 link modes, ∝ N.

The heavy mode carries no strain dependence: its rest energy is on-site. It moves only by hybridising into link modes that hop. The whole Hamiltonian is homogeneous of degree 1 in N, so Q_N = ΔE identically (44B).

**Vacuum.** Each species is at its 40A half filling, 3 of 7 modes per site, with the heavy band empty. m₀ sits inside the gap, so the heavy band is the lowest band above μ. Its detuning from the conduction bottom is Δ = e_c,min − m₀.

**Interaction.** On-site attraction between the heavy-mode densities of the two species, V = −U Σ_v n_h1,v n_h2,v. Its strain law is fixed now as w = 0: contact, independent of a, like the rest energy.

**Strain.** Q_a = Σ_links −a ∂H/∂a = t ∂E/∂t (isotropic). Directional versions use the x-links only (t_x ∂E/∂t_x), and likewise for y and z.

## Parts

### 45A — exciton total charge at its true minimum (Codex item 1; 44's model)

- Locate the dispersion minimum K* of 44's exciton at each U (E_b = 5 / 15 / 30% of gap), over the full K grid at L = 24.
- Then at K*, at L = 32:
  - E_b relative to the continuum minimum at that K;
  - total Q_a = L_bound + w⟨V⟩ for w = 0, plus the binding-only increment;
  - directional Q_a,x, Q_a,y, Q_a,z.
- Record the total-cancelling w = −L_bound/⟨V⟩ (descriptive).

### 45B — one heavy carrier

Grid: λ₀ = 0.1, m₀ ∈ {0.3, 0.5, 0.7} (gap −0.5 to 0.9, μ = 0.2), g ∈ {0.025, 0.05, 0.1, 0.2}.

For each point:
- the heavy band (Bloch, L = 32): bandwidth W_h, effective hopping t_eff from the band curvature at its minimum, and admixture 1 − |⟨h|u⟩|² at the minimum;
- at the band minimum: ΔE = e_min, Q_a = t ∂e_min/∂t (finite difference), and Q_a/ΔE;
- a moving state at k = (π/4, 0, 0): Q_a/ΔE, and its kinetic part (e(k) − e_min).

### 45C — bound pair of heavy carriers: binding stress

- The two-species two-body problem in the heavy band, at total momentum K = 0 and at its dispersion minimum. Separable on-site form, as in 44D, with species form factors.
- U is set so that E_b = 5 / 15 / 30% of W_h, at three (m₀, g) points spanning the 45B range.
- Measured:
  - ΔE_pair;
  - Q_a,pair;
  - the binding stress δQ_a = Q_a,pair − 2 Q_a,single(edge);
  - the lattice Nordtvedt ratio η = δQ_a / E_b;
  - the kinetic scale ⟨T⟩ = ⟨D⟩ − 2 e_min, so δQ_a/⟨T⟩.

## Predictions (scored)

- **P1 (45A).**
  - K* ≠ 0 (44 saw E₀ fall away from K = 0).
  - The total Q_a/ΔE of the exciton is > 1 at every U and at both λ₀: it stays radiation-like.
  - The total-cancelling w varies by > 20% across U (Codex's 1.9–3.7 at K = 0).
- **P2 (45B: heavy matter can be dust-like).**
  - Q_a/ΔE at the band minimum falls monotonically as g decreases at fixed m₀.
  - It scales as g² within 20%, comparing g = 0.025 with g = 0.05.
  - At g = 0.025 it is below 1e-2 at every m₀.
  - t_eff ∝ g² as well, so mobility and spatial charge vanish together.
  - The ratio R = Q_a(band min) / (bandwidth W_h) is recorded, not predicted.
- **P3 (45C: binding energy does not gravitate as energy).**
  - The lattice Nordtvedt ratio |η| = |δQ_a/E_b| exceeds 0.1 at every point. Binding carries a spatial charge of order its own size.
  - η is not constant: its spread across the 9 (E_b, (m₀, g)) points exceeds 20%.
  - Reason: with w = 0, δQ_a = δ⟨hop⟩ is set by how binding reshapes the pair's hopping. Nothing forces that to vanish, or to be a fixed fraction of E_b.
  - If |η| < 0.01 at all points, that is a surprise: binding would be stress-free on the lattice, the first lattice von Laue behaviour. The note would say so.
- **P4 (descriptive).** δQ_a/⟨T⟩, compared with the continuum virial expectation for a contact-bound pair. There is no clean continuum value for a contact pair, so this is recorded only.

## Checks

- **C1.** Homogeneity: Q_N = ΔE for one heavy carrier and one pair, within 1e-8 (uniform-N finite difference).
- **C2.** At g = 0 the heavy band is flat at m₀ and Q_a = 0 exactly (to 1e-12).
- **C3.** A dense projected-matrix solve at L = 10 equals the separable pair solver within 1e-10.
- **C4.** The interaction Hellmann–Feynman U ∂E₀/∂U equals ⟨V⟩ within 1e-8.
- **C5 (Codex's direct-derivative point).** For one pair, the hopping derivative of the two-body energy is computed with both particles allowed in all four bands above μ, at L = 4 (65,536 pair states, sparse Lanczos). It is compared with the heavy-band-projected value at the same L. The difference is reported as the projection error. Gate: below 5% of δQ_a.

## Outputs

- `src/exgc45.py` and `results/EX-GC45/`.
- A separate recomputation, `src/exgc45_recompute.py`, recomputes 45B Q_a/ΔE from eigenvalues of an independently built 7-band Bloch H, and η at one 45C point from a dense L = 10 solve.
