# EX-GC56 — the link shift as a substrate drift: one directed offset carries light and matter alike

2026-11-28 | EXPLORATORY | Claude | PI go-ahead 2026-11-28 01:25 UTC ("ok, go ahead with the substrate shift run"). Cap: 2 core-h (expected about 0.4). This plan is pushed before any run. Our construction; not LaViolette's Model G and not Chain Theory.

## Why

55 found that the operator which both sources and feels the shift has to be the lattice momentum P (53's law, GR's ADM form −β·p). A reciprocal energy current fails. But P was chosen because it is the conjugate of the coupling, not derived from the substrate. 56 asks what P means at the event level, and whether that meaning can be tested in a way that could fail.

## The idea (worked before any run)

**1. The shift coupling is literally transport of the pattern.**
- With P = Im(ψ_x†ψ_{x+a}) = (T − T†)/2i, where T is the lattice translation ψ(x) → ψ(x + a), the shift term in the evolution is −iH_βψ = (β/2)(Tψ − T†ψ). That is a real, spin-blind, mass-blind transfer of amplitude across each link, in one direction. It is classical advection, ψ_t = β∂ψ.
- So −β·P means that each directed link carries a fraction β dt of the pattern per tick. The substrate moves, and the pattern is carried with it.
- Light's shift term in the covariant 55 rule (u_t ⊃ β∂_su) is the same advection.

**2. One offset must carry everything equally (an equivalence principle for the shift).** If the shift is a drift of the substrate, every field must be carried by the same amount. On the lattice, both light and matter then shift their group velocity by Δv = −β cos k at wavenumber k, for every mass, both spins, and with or without the Wilson term. That is because P's spectrum is sin k, the same central-difference translation light uses.

**3. The energy current breaks this.** A coupling to J_E = ½{Ĥ, v̂} has spectrum E·v = ½∂_k(E²), so the drift is Δv = −β · ½∂_k²(E²).
- Naive fermion (r = 0): −β cos 2k. This is mass-independent, but not light's cos k.
- Wilson fermion (r = 1): −β[cos 2k + r² sin²k + (m + r(1 − cos k)) r cos k], which is 1 + mr at k = 0. It is mass-dependent.
- For light, ω·v_g = sin k, so for light the energy current coupling equals P. The two couplings differ only for matter.
- So J_E would make the substrate carry matter differently from light, and different masses differently. That is a lattice violation of the equivalence principle, and it is a cleaner reason than 55's well factor for why J_E fails.

**4. What is then derived and what is still put in.**
- **Fixed by the drift reading:** the coupling's normalisation. The coefficient of −β·P is 1 because matter must be carried exactly as light is.
- **From 54's counting argument:** the strength λ = 4 with which matter sources β.
- **Still an assumption:** the field law for β (coincidence transport like N's).

## Parts

**56A. Drift universality in uniform β** (β = 10⁻³ on all x-links, N = 1).
- **Matter.** Dirac packets with σ = 5, n = 48, T = 10, dt = 0.1, RK4. The centroid displacement with β, minus the same run with β = 0, over βT, against the packet average −⟨cos k_x⟩.
  - Cases: m = 0.3, 0.6 and 1.2; k₀ = 0, 0.15 and 0.3; r = 0 and 1; spin up and spin down (at m = 0.6, k₀ = 0.15).
  - Both couplings: P (53/55) and J_E (55's reciprocal operator, at N = 1).
- **Light.** The 55 covariant loop rule on an open stretch (L = 512), pulse width 16, k₀ = 0.25 and 0.5, T = 100. The centroid drift with and without β over βT, against −cos k₀ (with the packet-averaged value reported).

**56B. The conveyor identity** (cheap). For random ψ and random link β, check that −iH_βψ equals the real directed transfer (β/2)(Tψ − T†ψ), in its link form, to rounding. For J_E, check that the same quantity is not a pure translation: the residual against the best-fit advection is reported.

**56C. The Codex generalisation checks for the momentum rule.** The 55C pair setup (n = 49, T = 8, dt = 0.1, separation 14, P/P, λ = 4, all fields dynamic).
- **Cases:**
  - κ/2 at V = 0 and V = 0.2 (half the well depth);
  - V = 0.15 (with 55's V = 0 run);
  - packet width σ = 4 at V = 0 and V = 0.2.
- **Measure.** R_dyn = γF₁(V)/F₁(0), against 55's R_P = 0.9986 at κ, σ = 5 and V = 0.2.

## Predictions (scored)

- **P1 (matter drift, P coupling).**
  - Measured drift / (−β⟨cos k⟩) = 1 within 2 × 10⁻³ for every matter case: masses, k₀, r = 0 and 1, both spins.
  - The spread across cases is below 2 × 10⁻³. This is universal.
- **P2 (light drift).** Measured drift / (−β⟨cos k⟩) = 1 within 5 × 10⁻³ at k₀ = 0.25 and 0.5.
- **P3 (energy current is not universal).**
  - For J_E, the matter drift matches the band formula −β⟨½∂_k²E²⟩ within 1%.
  - At r = 1, m = 0.6 and k₀ = 0 it differs from P's drift by more than 30%, since 1 + mr = 1.6 against 1.
- **P4 (conveyor identity).**
  - For P, −iH_βψ equals (β/2)(Tψ − T†ψ) to 1e-14 relative.
  - For J_E, the best single-advection fit leaves a relative residual above 10⁻².
- **P5 (generalisation, momentum rule).**
  - |R_dyn − 1| < 1% in every 56C case.
  - |R(κ/2) − R(κ)| < 0.5% at V = 0.2.
  - |R(σ = 4) − R(σ = 5)| < 1% at V = 0.2.

## Checks

- **C1.** 56A with β = 0 in both runs gives zero drift to 1e-12.
- **C2.** 56A matter drift at dt/2 for one case (m = 0.6, k₀ = 0.15, r = 0, P) agrees within 1e-4 relative.
- **C3.** 56C at κ, σ = 5 and V = 0.2, recomputed through this study's code, reproduces 55's F₁ to 1e-10.

## What the outcomes mean

- **P1–P4 pass.**
  - The shift has a substrate reading: a directed drift of the medium, which carries light and every kind of matter with the same group-velocity shift.
  - That fixes the coupling to the momentum (lattice translation), with coefficient 1, by the equivalence principle rather than by choice.
  - The energy current is excluded on the same ground.
  - Together with 54's counting (λ = 4), the shift sector then has no free number. What is left by hand is the transport law of β itself.
- **P5 passes.** The momentum rule's Lorentz invariance of the co-moving pair holds beyond the single 53/55 calibration: half the well depth, another speed, another packet width.
- **P1 fails for some mass, spin or r.** The lattice momentum is not a universal translation for Dirac matter, and the drift reading breaks down.
- **P5 fails.** The pair result depends on the well or the packet, and 53's λ* has to be re-read.

## Outputs

`src/exgc56.py`, `src/exgc56_summary.py`, `src/exgc56_recompute.py`, `results/EX-GC56/`, note, LOG row and HANDOFF section.
