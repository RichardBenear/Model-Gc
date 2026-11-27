# EX-GC51 — one declared model end to end: Dirac matter, a coincidence-carried slowness field, and light

2026-11-27 | EXPLORATORY | Claude | **Plan only. Nothing is run until the PI gives a separate go-ahead.** This follows the PI's request of 2026-11-27 17:00 UTC ("go ahead with the third write up"), which is recommendation 3 of the [Codex review of 48/49/50](../review/EX-GC48_49_50_codex_review.md). Cap requested: 2 core-h. Our construction; not LaViolette's Model G, not Chain Theory.

## Why

Since 44 the field equations have been prescribed. 48, 49 and 50 each supplied one piece, under different assumptions: relativistic matter (48), a rule that carries N between neighbours (49), and a vacuum statement (50). Codex's point is that no single model has yet produced a source, its field, and a probe's motion from the same rules and parameters. 51 declares one such model, fixes every choice before any run, and tests it end to end.

## The declared model

Everything lives on one cubic lattice (h = 1). There is one scalar field N(x), the local firing rate. Each site advances in its own time τ, with dτ = N dt.

1. **Light** follows the 37B coincidence rule. A light field u obeys ∂_τ² u_x = Σ_y N_xN_y (u_y − u_x).
2. **Matter** is a single Wilson–Dirac fermion (48), treated first-quantized, with no Dirac sea.
   - Bond operators (kinetic hop plus Wilson bond, W1) carry the coincidence weight N_vN_w.
   - The mass term m β on a site carries N_v.
   - At uniform N this gives E² = N²m² + N⁴k² + …. The rest energy redshifts as N, and the limiting speed is N², the same as light.
   - **Lapse charge.** The charge density is q_x = ⟨∂H_m/∂N_x⟩, so Q = E + ⟨H_bond⟩. For a wide packet at rest under W1, ⟨H_bond⟩ → 0 (48B). So free dust has Q → E, which removes the order-1 equivalence-principle failure 43 found with the pair network.
3. **Field.** N is carried by the same coincidences. The transported quantity is the coordinate slowness of light, f(N) = −1/(2N²):
   (1/N)∂_t[(1/N)∂_t f(N)_x] = Σ_y N_xN_y [f(N_y) − f(N_x)] − κ q_x,  with κ = 4πG.
   - **Microscopic reading (declared, not derived).** Coordinate light delay per hop is 1/N². A coincidence passes on the delay it caused, not a rate, and delays add along a path.
   - 49 found that this choice gives β = 1. Codex notes it is not unique: β fixes only f″(1) = −3. 51 does not settle uniqueness, and it is listed as open.
4. **Back-reaction.** Matter moves in the evolving N, and N is sourced by matter's q. Nothing is prescribed except the parameters (m, r, κ, lattice size) and the initial states.

**What the model leaves out, by declaration.**
- The Dirac sea and its vacuum response. 50 showed that with coincidence weights the sea's response reappears in the clock sector, of order the vacuum link energy. 51 therefore uses first-quantized matter, and the cosmological-constant-type problem stays open.
- Bound composites. 48 showed contact binding fails, and they are deferred.
- A separate spatial factor a. Light's spatial half comes from the coincidence count, as in 37B.

## Parts

- **51A. Light speed of the field.**
  - **Setup.** At t = 0, a localized static fermion source (m = 0.5, σ = 3, frozen) is switched on in a flat N = 1 box, L = 96.
  - **Field measurement.** Measure the arrival time of half the final δN at distances r = 10–40 along the axis and the diagonal. The slope gives the field's signal speed c_N.
  - **Light comparison.** Do the same for a light pulse from a point source of u, with a matching spectral content.
- **51B. Static exterior.**
  - **Solve.** Solve the static field equation for a frozen fermion source at four strengths, with U(r = 12) between 0.01 and 0.05. Box: 129³, Dirichlet Gauss-law faces as in 49.
  - **Fit.** Use the registered 49 post hoc method, now prospective: fit the axis, face diagonal and body diagonal separately, and separate the source-linear part over strength.
  - **Outputs.** G_eff = A/(κQ/4π) and β.
- **51C. Probe fall (equivalence principle for free fermions).**
  - **Source and probes.** A static frozen source (from 51B) is placed with free probe packets released at rest at r = 16, with σ = 4. Probe masses m = 0.1 and 0.3; Wilson r = 1, 0.5 and 0 (naive fermion, k ≈ 0 packet, doublers unoccupied).
  - **Measurement.** Take the centroid acceleration from differential runs, with the source minus without, over t = 0–T. Compare with a_N = −∇N at the probe, taken from the same static solution.
  - **Prediction basis.** From 48, inertia is m/(1 + rm) and charge is m. So a/a_N = 1 + rm, and the lattice's only free-fermion EP violation is r Δm, of order mass/cutoff.
- **51D. Two free fermions, fully dynamic (reciprocity).**
  - **Setup.** Two packets (m = 0.3 and 0.5, σ = 4, separation 16, L = 64) and the field evolve together from rest. The field's initial state is the static solution of the pair.
  - **Measurement.** Measure the accelerations a₁ and a₂ after the retardation transient (t > separation/c). Check the momentum balance m_i,1 a₁ + m_i,2 a₂ ≈ 0, with m_i from 48's dispersion.
  - **Open point.** The coupled rule is not known to follow from an action, so a momentum imbalance is a real possible outcome.
- **51E. Light through the self-consistent field.**
  - **Setup.** Use the static 51B field at the largest strength. Pass light packets at impact parameters b = 8, 16 and 24.
  - **Measurement.** Measure the Shapiro delay relative to the one-site control rule (speed N) in the same field. 37B's ratio was 2 for a prescribed N; here N is the model's own field.

## Amendment A1 (2026-11-27, PI go-ahead 17:06 UTC, written with the code and before any run)

1. **51C comparators and probes.** A wave packet with internal momentum spread falls at a rate that depends on its velocity distribution (in GR too, since light falls twice as fast transversely). So 1 + rm holds only as σm → ∞. Three changes follow.
   - **Probes.** Use σ = 3/m with masses 0.25 and 0.5, so both have the same σm = 3. Evolve them in a 96³ sub-box of the 51B field at U(12) = 0.035, and run for T = 12.
   - **Comparators.**
     - a_GR = −N³⟨∂ₓN⟩, the coordinate acceleration of a slow particle in the model's metric g₀₀ = −N², g_ij = N⁻²δ_ij. This replaces the plain −∇N.
     - a_semi = −⟨∂²E/∂kₓ² · ∂E/∂N⟩⟨∂ₓN⟩ from the lattice dispersion, weighted over the packet's k distribution.
   - **Scored P3.**
     - a/a_semi = 1 ± 0.03 for every (m, r).
     - Within each r, a(0.25)/a(0.5) equals a_semi's ratio within 1%.
     - At r = 0, a(0.25)/a(0.5) = 1 within 1%.
     - Reported but not scored: a/a_GR, and a_semi/a_GR against 1 + rm.
2. **51D.**
   - **Initial state.** The field starts from the static solution of the frozen pair, so there is no retardation transient. Accelerations are fitted from t = 0.
   - **Reciprocity, scored P4.** The canonical lattice momentum ⟨kₓ⟩ of each packet gives the force F = dP/dt. The score is |F₁ + F₂|/|F₁| < 0.05.
   - **Diagnostic only.** The inertia-weighted acceleration balance.
   - **Parameters.** Masses 0.4 and 0.6, σ = 3/m, separation 20, n = 97, U(20) = 0.01.
3. **51A.** Use a flat Dirichlet boundary for both the dynamic run and its static target. The light comparison is the linear field u on N ≡ 1 with the same source.
4. **C4 scope, for the compute cap.** Halve dt for 51C (one case) and for 51D. Run the box 1.5 times larger for 51C only.

## Predictions (scored)

- **P1 (one cone).** c_N = 1.00 ± 0.02. Its arrival-time slope equals the light pulse's within 1%.
- **P2 (static).** G_eff/G = 1 within 1% at every strength. β = 1 ± 0.05 on each fit direction, after the source-linear separation.
- **P3 (fall).** a/a_N = 1 + rm within 3% for each (m, r). Within each r, the ratio of the two masses' accelerations equals (1 + r m₁)/(1 + r m₂) within 1%. At r = 0 the two masses fall equally within 1%.
- **P4 (reciprocity).** |m_i,1 a₁ + m_i,2 a₂| / (m_i,1 |a₁|) < 0.05. This is diagnostic. If it fails, the combined rule does not conserve momentum, and the next step is to find the action, or show that none exists.
- **P5 (light).** The Shapiro ratio to the one-site control is 2.00 ± 0.06 at each b.

## Checks (separate paths)

- **C1.** The 51B static solver with a Gaussian source reproduces 49's slowness solution to 1e-10.
- **C2.** The matter Hamiltonian at uniform N = 0.9 and 1.1 reproduces the scaled 48 dispersion, E(k; N) = N·E_48(k) with kinetic ∝ N², to 1e-12.
- **C3.** The charge identity Q = E + ⟨H_bond⟩ holds by direct finite difference in uniform N, to 1e-9.
- **C4.** For 51C and 51D, report the largest energy excursion during the run, not only the end drift, and repeat with dt halved and with L larger by 1.5. The fitted accelerations must agree within 1%.

## What the outcomes mean

- **P1–P3 and P5 pass.** Weak-field gravity (Newton's law, light bending with γ = 1, β = 1, and the equivalence principle for free fermions up to an r m lattice term) comes out of one set of rules, with no field equation put in. That would be the programme's first end-to-end mechanism candidate. Three questions would remain: why slowness, the vacuum, and composites.
- **P4 fails.** The model is not yet a consistent dynamics, because action and reaction are unequal. That points to deriving the field rule from an action that includes matter's coupling.
- **P3 is off by more than r m.** Something besides the Wilson term breaks universality, and I will trace it.

## Outputs

`src/exgc51.py`, `src/exgc51_recompute.py`, `results/EX-GC51/`, note, LOG row and HANDOFF section.
