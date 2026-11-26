# EX-GC37A — does gravity's stiffness come from the matter field? (induced gravity on Y1)

2026-11-25 | EXPLORATORY | Claude | PI go-ahead in thread ("go ahead with both", 2026-11-26 00:13 UTC). Written before any run. Default cap 2 core-hours; expected well under 0.5. This is our construction, not LaViolette's Model G.

## Why

Y1 reproduces linearised GR (34N–35K) only because its vacuum action *is* linearised Einstein–Hilbert (EH), put in by hand. That gives no mechanism. Sakharov's induced-gravity idea (1967) supplies one: the stiffness of spacetime, Newton's 1/G, is the vacuum response of the matter fields that live on it. Space resists bending because bending it changes the zero-point energy of everything in it. If that holds on Y1, "mass curves spacetime" and "spacetime resists curving" have the same origin: the matter substrate.

The question is whether this actually happens on our lattice, with the right structure. In the continuum the answer is yes by symmetry. A Lorentz-invariant cutoff can only induce EH at order q², with a cosmological term, a graviton mass to tune away, and a positive 1/G for a minimally coupled scalar. Y1's cutoff is not Lorentz invariant: space is a lattice, time is continuous, and the matter dispersion is ω = |K| with K_a = 2 sin(q_a/2)/h. The induced q² term is dominated by modes near the cutoff. So the lattice may induce something gravity-like but not EH, for example with a graviton speed different from light or an extra scalar mode. That is the test.

## Question

Integrate out one lattice scalar field on Y1 in the vacuum, at one loop. Does the induced quadratic action for the metric perturbation (lapse n, shift β, h_ab) have:
1. the correct sign (attractive, positive-energy gravitons);
2. the linearised-EH structure at order q²;
3. so graviton speed = light speed, exactly two propagating modes, and static γ = 1?

## Method

**Matter action to second order.**
- 34N/34O's `mlagrangian` couples φ to the metric only at linear order. Its variation gives the stress tensor, but the induced action needs the contact (seagull) terms too.
- So write one lattice matter action as an explicit function of the metric: √γ/N φ̇² minus N√γ γ^{ab} G_a G_b, with the shift term, on Y1 placements.
- Expand it to second order in (n, β, h).
- It must reduce to 34O at first order, including the product-form normal stress.
- Fix the second-order placement choices once and write them in the code before any spectral run. No tuning after results.

**Induced kernel.** Use the standard one-loop formula in Hamiltonian (spectral) form. Space is a lattice and time is continuous, so every sum is finite.
- Π(Ω, q) = ⟨T T⟩ two-particle spectral sum + ⟨δ²L/δh δh⟩ contact term.
- Evaluate it on N³ Brillouin-zone grids, N = 32, 64, 128, and confirm convergence.
- Work in Euclidean frequency Ω and continue analytically for the wave readouts.

**Vacuum subtraction.** The zeroth-order term (vacuum energy, which acts as a cosmological constant and graviton mass) is reported, then subtracted as one counterterm. That is the one standard tuning in Sakharov's scheme, and it is recorded as such.

**Readouts.**
- Small-(Ω, q) expansion of the subtracted kernel at q², split into the lapse, shift, trace and transverse-traceless (TT) sectors.
- Compare with the 34N EH kernel, sector by sector, and fit c_g² (graviton speed²), λ_ind (the trace kinetic coefficient; EH has λ = 1), and the static lapse/space ratio γ.
- Sign and size of 1/G_ind in lattice units.
- Repeat at scalar mass m = 0, 0.5/h and 1/h. Induced 1/G should be set mainly by the cutoff and fall slowly with m.

**Stages.**
- **S1, TT sector only.** This is gauge-safe at linear order and needs only the gradient-term contact piece. It gives the sign of G_ind and c_g.
- **S2, full sector.** Gives λ_ind, γ and the mode count.

## Internal checks (gates)

- **G1, Ward identities.** The induced kernel must be annihilated by the Y1 re-slicing and displacement gauge directions, to roundoff, at the four orders above zero, after subtraction.
  - If the second-order action cannot be made exactly gauge invariant (the 34S–34Z experience), G1 is scored as "residual / |Π₂| → 0 as (qh)²".
  - S2 is then reported with that qualifier.
- **G2, continuum control.** The same code, with the lattice symbols replaced by continuum momenta and a smooth Lorentz-invariant cutoff, must reproduce:
  - pure EH structure at q²;
  - c_g = 1, γ = 1;
  - positive 1/G.
- **G3, separate path.** Independent finite-difference evaluation: second derivative of the vacuum energy of the lattice Hamiltonian in a static, weak, Fourier-mode metric background. Gives the static sector at small q and must agree with S2's static kernel.

## Predictions (scored)

- **P1:** 1/G_ind > 0 at m = 0 (attractive gravity). Expectation: yes.
- **P2:** TT sector at q² ∝ (Ω² − c_g² q²), with |c_g² − 1| ≤ 2%. Expectation, recorded and not scored: **fails**, with an O(10%) lattice deviation, because the induced term is dominated by cutoff modes.
- **P3:** static γ_ind within 2% of 1. Expectation: uncertain.
- **P4:** exactly two propagating modes, so λ_ind = 1 within 2%. Expectation: fails.

## What each outcome means

- **All pass:** the matter substrate alone produces GR's linear stiffness on Y1. This is a candidate mechanism: gravity as the vacuum elasticity of the discrete matter field. The next step is a dynamical run with the induced action replacing the hand-written one.
- **Right sign, wrong structure (expected):** the mechanism exists in principle, but the lattice cutoff breaks the GR structure. That points to what the substrate needs: a cutoff that respects local Lorentz symmetry, i.e. discrete time too, or a random or causal-set-like lattice. That is a concrete design criterion, not a dead end.
- **Wrong sign:** pure scalar induced gravity is ruled out on Y1. Fermion or gauge-field content would be needed, as in the continuum.

## Outputs

`src/exgc37a.py`, `src/exgc37a_recompute.py` (G3), `results/EX-GC37A/` (note and JSON).

Before running: an `explore/LOG.md` entry and PI go-ahead.
