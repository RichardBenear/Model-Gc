# EX-GC47 — two bound heavy pairs: do they attract through the prescribed fields, and how do they fall?

2026-11-27 | EXPLORATORY | Claude | Plan only. PI request in thread: "the review of 46 is available. I would like to see the tests outlined there done." (2026-11-27 04:35 UTC). That message is taken as the run go-ahead for the tests the [Codex review of GC46](../review/EX-GC46AB_codex_review.md) outlines (its "Proposed GC47", stages 1–3). The plan is pushed before any run. Cap 2 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why

46 produced a bound pair of heavy carriers with small spatial charge (Q_a/E ≈ 0.8% at g = 0.025) and a small but finite mobility. Codex proposes using two such pairs, one the source and one the probe, with no direct matter interaction between them, and asking whether the prescribed fields make them attract: first as a static cross energy, then as motion.

What this can and cannot show (stated before any run). The field equations, the vacuum prescription and the way the fields couple to matter are all put in. They are the 44 equations, not derived from the network. The test is whether the *microscopic* matter states, with their own charge profiles and inertia, give a reciprocal attraction with the expected coefficient, and how the bound bodies then move. It does not derive gravity from the network, and it does not test universality of free fall beyond the two benchmarks.

**A prediction made before running, which the plan is built to check.** The pair's inertia is set by its hopping curvature, m_i = 1/(2 t_pair). The 46 values (t_pair ≈ 6e-5 at the first benchmark, pair energy E ≈ 0.99) give m_i ≈ 8,500 against an energy of about 1. The gravitational charge is the energy, not m_i. So the probe's acceleration should be 2 t_pair times the Newtonian value, about 10⁻⁴ of it, and it should differ between pairs with different t_pair. That would be a free-fall failure of a new kind: not from the charge, but from inertia that is not tied to energy. Making Q_a/ΔE small (46B) makes t small at the same rate, so dust-likeness and a Newtonian acceleration cannot be had together on this lattice.

## Model

- **Matter.** The 45/46 heavy-carrier pair: two species, on-site contact V = −U Σ n_h1 n_h2, lapse-homogeneous weights, hop ∝ 1/a. Carriers in the heavy band (index 3 of 7) only for the charge profiles and dynamics. 46 showed the four-band correction to E_b is 0.3% at weak hybridisation.
- **Benchmarks.** B1: (m₀, g) = (0.5, 0.025), E_b/W_h = 0.15. B2: (0.7, 0.025), E_b/W_h = 0.15. λ₀ = 0.1.
- **Fields (44, prescribed).** Lapse n_v = N_v − 1 = Φ_v on sites. Spatial factor a_vw = exp(−(Ψ_v + Ψ_w)/2) on links. At linear order, ∇²Ψ = 4πG ρ_N and ∇²Φ = 4πG (ρ_N + ρ_a), with the local vacuum counterterm (43's prescription L), which divides both by 1 + 4πG c₂χ, c₂χ = 1.001 (43). 4πG = 0.01 (main), 0.005 and 0.02 (linearity). Static source, periodic box, zero mode removed.
- **Densities.** ρ_N,v = ⟨H_site,v + ½ Σ_links at v H_link + V_v⟩, with V_v the contact term at v (it carries the lapse like any site term). ρ_a,v = ½ Σ_links at v ⟨H_link⟩. Σρ_N = E (homogeneity) and Σρ_a = Q_a.
- **Probe coupling.** δE_P = Σ_v [Φ_v ρ_N,P(v) + Ψ_v ρ_a,P(v)]. This is first order in the probe's coupling to the source-generated fields, and it counts the interaction once. Bare far-field prediction: V₁₂ = −G (E₁E₂ + E₁Q_a2 + Q_a1E₂)/r.
- **Pair form factors.** The pair at total momentum K = 0 has relative amplitude φ(k) ∝ Φ(k)/(D(k) − E). Densities in momentum space are the transition matrix elements F_N(q), F_a(q) between the pair states at K and K + q, built from single-carrier matrix elements of the lapse and strain operators between Bloch states at k and k + q (bond weight e^{ik·d}(1 + e^{iq·d})/2), plus the contact part. The dependence of φ on K is neglected: the COM energy scale (t_pair K²) is far below E_b. A packet with Gaussian envelope of std σ has ρ(q) = F(q) exp(−σ²q²/2) exp(−iq·R). F(q) depends only on the sorted |n| of q by cubic symmetry, so it is computed on the wedge.

## Parts

### 47A — benchmark pairs (stage 1 of the review)

At L = 32, for B1 and B2:
- E(K) along x, y, z and (1,1,1) at K = 2πn/32, n = 0..8. Fit E(K) − E(0) = t K² + u K⁴ to get t_pair; isotropy; comparison with 46's one-step value.
- Pair size: rms relative displacement from the relative wavefunction. Gap to breakup (E_b, and its K dependence).
- E, Q_a, E_b; Q_a/E; the inertia ratio m_i/E = 1/(2 t_pair E).

### 47B — static cross energy (stage 2)

- F_N(q), F_a(q) for B1 and B2 on the wedge at L = 16, 24 and 32. Checks F_N(0) = E, F_a(0) = Q_a.
- V₁₂(R) for R along x, y, z at separations 3–12 (L = 32; the smaller L at shorter range), for packets σ = 0 and 2, source B1 / probe B1, and B1/B2 both ways.
- Comparators: point charges (E, Q_a) in the same periodic Green function; the bare dust value (E² only).
- Controls: G = 0; 4πG = 0.005/0.01/0.02 (linearity); a one-body control (one carrier as source and probe); swap source and probe.

### 47C — release (stage 3)

- **Reduction.** The probe's centre-of-mass wavefunction c(R) evolves under H_eff = E(K̂) + V₁₂(R̂ − R_S), with E(K) the exact pair band and V₁₂ the σ = 0 cross energy of 47B. This is first-order degenerate perturbation theory in the fields on the bound branch. It is valid while the field varies by much less than the gap E_b across the pair, and while second-order terms are small. Both are checked.
- **Validation of the reduction.** At L = 6, the exact two-body problem (both carriers in the heavy band, full pair space of 46,656 states, sparse Lanczos) in the same static source field is compared with H_eff on the 216 COM states: lowest 5 levels.
- **Probe release.** Source fixed, probe packet σ = 3 at rest, separations 6 and 8 (L = 24, 3D FFT split-step, dt convergence by halving). A source-off run is subtracted. Measure the initial acceleration of ⟨R⟩ and compare with −2 t_pair ∂V₁₂/∂r and with the Newtonian value G_eff E_S/r².
- **Both released.** Relative coordinate at total momentum zero, H_rel = E₁(k) + E₂(−k) + V₁₂(r), for B1/B1 and B1/B2. Report the approach, and reciprocity.

## Predictions (scored)

- **P1 (47A).**
  - t_pair from the K → 0 fit is within 10% of 46's single-step value at both benchmarks.
  - E(K) is cubic-isotropic: the spread of t over x, y, z, (1,1,1) is below 5%.
  - The pair is compact: rms relative separation below 5 sites.
  - m_i/E exceeds 10³ at both benchmarks.
- **P2 (47B).**
  - V₁₂ < 0 at every separation, packet width and ordering.
  - At r ≥ 6, V₁₂ is within 3% of the point-charge periodic prediction, for σ = 0 and σ = 2.
  - V₁₂ / V_dust is within 3% of 1 + 2Q_a/E (about 1.0166) at r ≥ 6 for B1/B1.
  - Reciprocity: V₁₂(B1 source, B2 probe) = V₁₂(B2 source, B1 probe) to 1e-10 relative.
  - V₁₂/(G_eff) is constant to 1e-9 between the three 4πG values, and V₁₂ = 0 at G = 0.
  - F_N(0) = E to 1e-8 and F_a(0) = Q_a to 1e-6.
- **P3 (47C).**
  - H_eff reproduces the exact two-body lowest 5 level shifts to within 3%.
  - The measured initial acceleration agrees with −2 t_pair ∂V₁₂/∂r within 5%.
  - a/a_Newton is below 10⁻³ at both benchmarks (it is about 2 t_pair E_P(1 + …)/E_P).
  - a(B1)/a(B2) differs from 1 by more than 20%: free fall is not universal between these two bound bodies.
  - Both released: the separation decreases for B1/B1 and B1/B2; dt halving changes the displacement by less than 1e-3 relative.
- **P4 (descriptive).** The box size at which V₁₂ r/G stops changing by more than 1% at fixed r. The deviation from the point-charge value at r < 6, which shows the pair-size (form-factor) correction.

If a part is infeasible within the cap, the note will say which stage was not done.

## Checks

- **C1.** The form-factor sums at q = 0 reproduce E, ⟨V⟩ and Q_a of the solver (above).
- **C2.** The wedge F(q) equals a direct sum over a random sample of full-grid q (cubic symmetry used correctly), to 1e-10.
- **C3.** Exact vs H_eff at L = 6 (P3 first clause).
- **C4.** Norm and energy conservation of the COM evolution to 1e-10.

## Outputs

- `src/exgc47.py` and `results/EX-GC47/`.
- A separate recomputation, `src/exgc47_recompute.py`: direct real-space construction of the pair density operators at L = 6 (own 7-band build), cross energy computed by summing real-space densities against the real-space field, compared with the form-factor result.
- Not done here: deriving the field equations from the network, free-fall universality over more than two bound bodies, a finite-range binding interaction (Codex GC45 test 3, still deferred), and the frozen valence sea.
