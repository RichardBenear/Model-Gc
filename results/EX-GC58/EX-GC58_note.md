# EX-GC58 — the shift as local events: directed transfers on links, circulation on plaquettes, and the spin term GR needs

2026-11-28 | EXPLORATORY | Claude | PI closure pending | [Plan](../../spec/EX-GC58_plan.md) (pushed before any run) | About 0.62 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3; 4 workers). Scoring: `src/exgc58_summary.py` → `summary.json`; recomputation: `src/exgc58_recompute.py` → `recompute.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**The link shift's coupling can be written as two kinds of local event, and with both of them matter feels and sources frame dragging as GR says, spin included.** All six registered predictions and all four checks passed.

- **Link event:** a directed coincidence on a link moves amplitude one way. It rotates the two sites' amplitudes by β dt/2, the same for every spin and mass.
- **Plaquette event:** the circulation of the offsets round a plaquette rotates the spinors at its corners by a quarter of the local curl (dt/4 times the curl).
- **The source:** what each link records is the conjugate count of what its events did to matter, P + ¼Cᵀs.

What each part found:

- **It is a rule (P1 PASS).**
  - Each event preserves the norm to 7 × 10⁻¹⁶.
  - A sweep of events matches the exact evolution, and four random event orders match each other, up to errors that scale exactly as (β dt)². The fitted slopes are 2.0000 for links alone and for links with plaquettes.

  So the order of events does not matter at first order: the rule is causally invariant to O(β² dt²) per tick.
- **Drift universality survives the event form (P2 PASS).** Across ten cases (masses 0.3–1.2, r = 0 and 1, k₀ = 0–0.3, both spins), event stepping carries matter by −β⟨cos k⟩ to 7 × 10⁻⁵, and agrees with 56's Hamiltonian runs to 2.2 × 10⁻⁵.
- **Spin now feels GR's frame dragging (P3 PASS).** In 57's rotating-ring field, a resting spin packet's precession per unit of total angular momentum (ΔE/ΔJ_z) is:

  | | σ = 3 | σ = 5 |
  |---|---|---|
  | Links and plaquettes | 1.00014 | 0.99971 |
  | Links only | 0.099 | 0.046 |

  These are fractions of the local Ω = −½∇×β; d = 32 gives the same to 2 × 10⁻⁴. The spin-balanced orbiting probe changes by at most 1.1 × 10⁻⁴, and in a pure gradient the plaquette term is zero to rounding.
- **Spin gravitates (P4 PASS).** A resting spin-up packet with links and plaquettes sources a shift whose angular momentum is L + S = 0.4988, so its j = ½ is fully counted, to 2 × 10⁻¹³. With links only, the source is L = 0.049. Far loops read the dipole of that J to 0.1%.
- **Still one interaction (P5 PASS).** With the spin term, the energy drift is 3.3 × 10⁻⁵ of the interaction energy and falls 32× at dt/2 (RK4).
- **Lorentz invariance kept (P6 PASS).** The co-moving pair with the spin term gives R = 0.9967, against 0.9986 without it. Fitting t ≥ 2 only gives 0.9985. The spin terms change F₁ by −0.04% at V = 0 and −0.24% at V = 0.2.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | norms < 1e-13; error and order spread slopes 2.0 ± 0.15 (P and PS) | 6.7e-16; 2.0000, 2.0000 | PASS |
| P2 | event drift / (−β⟨cos k⟩) within 2e-3; vs Hamiltonian within 1e-4 | 7.0e-5; 2.2e-5 | PASS |
| P3 | PS ΔE/ΔJ / Ω within 3%; P < 0.15; orbit change < 1%; gradient < 1e-12 | 0.9997–1.0001; 0.046–0.099; 1.1e-4; 3e-36 | PASS |
| P4 | J_PS = L + S, J_P = L within 1%; loops / dipole within 5% | 1.0000, 1.0000; 0.9993–1.0009 | PASS |
| P5 | PS energy drift < 1e-4, falls ≥ 10× at dt/2 | 3.3e-5; 32× | PASS |
| P6 | \|R_PS − 1\| < 1.5%; \|R_PS − 0.9986\| < 0.5% | 0.9967; 0.19% | PASS |

## Checks

| | | |
|---|---|---|
| C1 | β = 0 events are the identity | 0 — PASS |
| C2 | H_pl Hermitian; ⟨H_pl⟩ = −Σβ·¼Cᵀs | 7.8e-15; 1.9e-15 — PASS |
| C3 | 58C links-only reproduces 57's spin ratios | exact — PASS |
| C4 | 58F links-only reproduces 55's F₁ at V = 0.2 | 0 — PASS |

**Recompute.**
- The vectorised event sweep, a different code path from the sequential event loop, gives slope 2.0000.
- Spin precession recomputed as an operator expectation, with the curl built on the probe box, agrees to 10⁻¹³.
- The spin part of the source's angular momentum equals ½⟨Σ_z⟩ to 2 × 10⁻¹³.
- The energy drift recomputed from the stored histories agrees.
- R from an independent least-squares fit is 0.99669, or 0.99847 for t ≥ 2.
- Shared dependencies: 53's ring, probe and shift code, and 55's evolve.

## Deviations

None in the scored runs. 58F's C4 run was not skipped: the whole study stayed at 0.62 core-h, under the 1.5 core-h threshold in the plan.

## What it means

- **The shift sector is now a local event rule with GR's coupling.**
  - A directed coincidence moves amplitude across its link.
  - A circulating set of offsets rotates spin.
  - The field records the conjugate counts of both.
  - Together these give Ω·(L + S): spin and orbit feel the same frame dragging, and both source it.
- **What the rule buys.**
  - Matter's coupling and source come from one generator, so energy balance needs no extra choice.
  - The co-moving pair stays Lorentz-invariant to 0.3%.
  - Event order does not matter at first order.
- **What is still put in:**
  - the strength 4 (54's two-site Doppler-squared counting);
  - the field law for β (how offsets propagate between links);
  - the plaquette event's weight (¼, the Belinfante value). It is fixed by requiring a resting spin to carry total angular momentum ½, which it then does to 2 × 10⁻¹³, but it is not derived from counting.
- **Suggested next (not approved):**
  1. A transport rule for the offsets themselves, so that β's wave equation comes from coincidences as N's does. Score it on 53's pair R and the ring's far field.
  2. Gravitational radiation from the full scalar-plus-shift system (deferred since 51).
  3. A counting argument for the ¼, for example from how many plaquettes share a site and how many links each plaquette rotates.
