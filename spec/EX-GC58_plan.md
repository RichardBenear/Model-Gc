# EX-GC58 — the shift as local events: directed transfers on links, circulation on plaquettes, and the spin term GR needs

2026-11-28 | EXPLORATORY | Claude | PI go-ahead 2026-11-28 ("Your 57 plan and 58 plan are very good. go ahead with them."). Cap: 2 core-h (expected about 0.5). This plan is pushed before any run. Our construction; not LaViolette's Model G and not Chain Theory.

## Why

- 55–57 fixed how matter feels and sources the link shift: the coupling is −β·P, the source is P with strength 4, and in uniform and rotating shifts this carries light and orbiting matter exactly as GR does.
- Two things are still open:
  1. the coupling is written as a Hamiltonian, not as local events, so it has not been tested as a rule (unitary per event, independent of event order);
  2. 57 found that an elementary spin feels only its orbital share of frame dragging (0.11 at σ = 3, 0.05 at σ = 5), because P is the canonical momentum density. GR couples g₀ᵢ to the symmetric (Belinfante) density, P_B = P + ¼∇×s with s = ψ†Σψ.
- 58 writes the shift sector as two local event types and tests whether they give the 55–57 physics plus the missing spin term.

## The rule (worked before any run)

1. **Link event (directed transfer).** A coincidence on link (x, x + a) with timing offset β moves amplitude one way. It rotates the pair (ψ_x, ψ_{x+a}) by the real angle θ = β dt/2, the same for every spinor component:
   - ψ_x → cos θ ψ_x + sin θ ψ_{x+a};
   - ψ_{x+a} → −sin θ ψ_x + cos θ ψ_{x+a}.

   This is exactly exp(−iH_β,link dt) for that link (56's identity), so each event is unitary on its own. A sweep over all links equals exp(−iH_β dt) up to commutators of order (β dt)². The order of events therefore does not matter at first order: the rule is causally invariant to O(β² dt²) per tick.
2. **Plaquette event (circulation).** The circulation of β round a plaquette is a local rotation of the substrate. It rotates the spinors at the plaquette's corners by exp(i(dt/4)(C β)·Σ), where Cβ is the site curl from 57 (the mean of the 4 plaquettes touching the site). As a Hamiltonian it is H_pl = −¼ Σ_x s(x)·(Cβ)(x), the Belinfante term. Then the interaction is H_int = Ω·(L + S), with Ω = −½∇×β, which is GR's spin–rotation coupling. It vanishes identically for β = ∇χ.
3. **Counting (the source).** What a link's offset records is the conjugate of what its events did to matter: the transfer count P_l from link events, plus ¼(Cᵀs)_l from the plaquettes the link borders. The source is therefore j = ∂⟨H_int⟩/∂β = P + ¼Cᵀs, and energy and momentum balance follow from the one generator. The strength 4 is still 54's counting, and the field law for β (the transport of offsets) is still put in.
4. **What spin should then do.** For a spin packet, ΔE/ΔJ_z = Ω, with J_z = L_z + ½⟨Σ_z⟩. Under links only, it is ΔL/ΔJ ≈ 0.1. A resting spin then sources a shift with angular momentum J = L + S, where under links only it sources L. A spin-balanced orbiting probe is essentially unchanged.

## Parts

**58A. The events as a rule** (n = 13, random ψ, random β of scale s·10⁻³, dt = 0.1, s = 1, 2, 4).
- Each link rotation and each plaquette spin rotation preserves the norm.
- One sequential sweep of link events (lexicographic order) is compared with exp(−iH_β dt)ψ (Taylor series to convergence). The error is reported against s.
- Four random event orders are compared with each other. Their spread is reported against s.
- Repeat with link and plaquette events together, against exp(−i(H_β + H_pl)dt).

**58B. Drift universality by events.** As in 56A, with β = 10⁻³ uniform on x-links, n = 48, T = 10, dt = 0.1.
- Matter is stepped by Strang splitting: half a sweep of link events (the six parity classes forward), then RK4 for H_m, then half a sweep in reverse.
- Cases:
  - m = 0.3 and 1.2, r = 0 and 1, k₀ = 0 and 0.3;
  - spin up and spin down at m = 0.6, r = 0, k₀ = 0.15.
- The drift is compared with −β⟨cos k⟩ and with the Hamiltonian (RK4, 56 code) drift of the same case.

**58C. The spin term in the 57 ring field** (53's n = 129 ring, ℓ = ±2, momentum-rule shift, static).
- Resting spin packets (σ = 3 and 5) on axis at d = 24 and 32: ΔE/ΔJ_z under links only (P) and under links plus plaquettes (PS), against the 57 curl Ω.
- 57's spin-balanced orbiting probe under P and PS.
- 57's gradient field β = ∇χ: the PS spin energy.

**58D. Spin sources frame dragging** (n = 97, a resting spin-up packet, σ = 3, at the centre).
- Source j under P and PS. Its angular momentum angmom(j) is compared with L_can and with L_can + ½⟨Σ_z⟩.
- The static shift is solved (53's `poisson`, dipole boundary from angmom(j)). Square loops of half-side 12 on axis at d = 20 and 28 give the loop readout over the dipole prediction for that J.

**58E. Energy bookkeeping with the spin term** (55B setup: N frozen, n = 49, one spin-up packet at V = 0.2, T = 8, dt = 0.1 and 0.05, dynamic β). Rule PS, with source P + ¼Cᵀs.

**58F. Co-moving pair with the spin term** (55C setup: n = 49, T = 8, separation 14, V = 0 and 0.2, all fields dynamic, both spins up along z). R_PS = γF₁(0.2)/F₁(0), against 55's R_P = 0.9986. The change in F₁(0) between PS and P (the static spin–spin part) is reported.

## Predictions (scored)

- **P1 (a local, order-free rule).**
  - Norms are preserved to 10⁻¹³.
  - The sweep-vs-exponential error and the order spread both scale as (β dt)²: fitted log–log slope 2.0 ± 0.15 over s = 1, 2, 4, for links alone and for links plus plaquettes.
- **P2 (drift universality survives the event form).** In every 58B case:
  - drift / (−β⟨cos k⟩) = 1 within 2 × 10⁻³;
  - agrees with the Hamiltonian drift within 10⁻⁴ relative.
- **P3 (the spin term gives GR's spin drag).**
  - Under PS, ΔE/ΔJ_z / Ω_curl = 1 within 3% at σ = 3 and 5, d = 24 and 32. Under P it stays below 0.15.
  - The orbiting probe's Ω changes by less than 1% between P and PS.
  - In the gradient field, the PS spin energy is below 10⁻¹² of the ring value.
- **P4 (spin gravitates).**
  - angmom(j_PS) = L_can + ½⟨Σ_z⟩ within 1%, and angmom(j_P) = L_can within 1%.
  - Loop readout / dipole(J = angmom(j)) within 5% at both distances, for both rules.
- **P5 (still one interaction).** Under PS, the energy drift over the interaction energy is below 10⁻⁴ and falls at least 10× at dt/2.
- **P6 (Lorentz invariance kept).** |R_PS − 1| < 1.5% and |R_PS − 0.9986| < 0.5%.

## Checks

- **C1.** With β = 0, all events are the identity (to rounding).
- **C2.** H_pl is Hermitian, and ⟨H_pl⟩ = −Σ_l β_l·¼(Cᵀs)_l to 10⁻¹² (source is the exact conjugate).
- **C3.** 58C under P reproduces 57's spin ratios (0.10943 and 0.04776 at d = 24) to 10⁻⁸.
- **C4.** 58F under P reproduces 55's F₁ at V = 0.2 to 10⁻¹⁰. This needs one extra run, and is skipped if the projected total exceeds 1.5 core-h.

## What the outcomes mean

- **P1 and P2 pass.** The shift's coupling is a local, unitary, order-free event rule: a directed coincidence moves a fraction β dt/2 of amplitude across its link, spin- and mass-blind. This is the event form 56 inferred.
- **P3–P6 pass.** One more local event, circulation round a plaquette rotating the spinors, completes the coupling to GR's symmetric momentum density.
  - Spin and orbit then feel the same frame dragging.
  - Spin sources it.
  - Energy balance and the co-moving pair's Lorentz invariance survive.

  What is still put in is the strength 4 (counting) and the field law for β.
- **P3 passes but P6 fails.** The spin term is GR-like for a probe in a given field, but it disturbs the pair. That would need a look at spin–orbit terms the momentum rule does not carry.
- **P1 fails.** The order of link events matters at first order. That would mean the shift cannot be read as independent local coincidences.

## Outputs

`src/exgc58.py`, `src/exgc58_summary.py`, `src/exgc58_recompute.py`, `results/EX-GC58/`, note, LOG row and HANDOFF section.
