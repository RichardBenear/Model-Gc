# EX-GC63 — the directed signal count, per local tick: does it supply the shift with nothing tuned?

2026-11-28 | EXPLORATORY | Claude | Plan only; not run. PI chose "Directed signal-count rule (3D)" on the decision card; this plan is the refined version proposed after finding that the chosen rule is 54's D2 rule (EX-GC62 note, corrected section 1b). Go-ahead for this refined version pending. Requested cap: 2 core-h (expected about 1). Our construction; not LaViolette's Model G and not Chain Theory.

## Why

EX-GC62 (section 1b) argues that λ = 2(1 + γ) is the no-preferred-frame condition, and that the shift should be the forward-minus-backward asymmetry of the signal count, with each endpoint Doppler-shifted: (1 ∓ v)² per direction, so 2 × (1 + γ) = 4 in 3D. That is 54's D2 rule. 54 measured it with the lattice energy current ½{Ĥ, v̂} as "count × velocity" and found:
- strength 3.85 in uniform N, against 53's momentum law at 4;
- strength 3.37 in place, inside the packets' own clock-rate wells, because Ĥ and v̂ carry hop weights N_x N_{x+a};
- R_dyn(0.2) = 1.019, against 1.002 ± 0.008 registered.

54 named the open question: count per local tick. A coincidence on link (x, x+a) occurs at coordinate rate N_x N_{x+a}; counted in the endpoints' own ticks, the weight is 1. So the per-local-tick count is the energy current with flat weights.

What this tests, and what it does not:
- **It tests** whether removing the clock-well factor makes the counted source keep the co-moving pair Lorentz-invariant with nothing tuned, and whether the count also carries a resting spin's source.
- **It does not** derive the two-site weight 1 + γ = 2 (that is 37B's, and the D study's question) or the Doppler square (argued in EX-GC62 1b, not derived from a rule here).

## The rule

Source on link (x, a): s = 4κ · ½ Re[ψ_x†(Ô ψ)_x + ψ_{x+a}†(Ô ψ)_{x+a}], with Ô = ½{Ĥ₁, v̂₁ₐ}. Ĥ₁ and v̂₁ are 54's `Hop` and `vop` evaluated with N ≡ 1 (mass m per local tick, hop weight 1). Code: `exgc54.rule_current(psi, ones, m, "D2E")`. Everything else is 53E/54B exactly (n = 65, separation 14, σ = 5, spin up, T = 16, dt = 0.1, source refreshed per step and held across RK4 stages as in 54's A1, initial shift from the boosted static solve). Rule label: **D2-E1**.

## Parts

- **63A (co-moving pair).** D2-E1 at V = 0 and V = 0.2, pair and single packet each; F₁ and R_dyn = F₁(V)/F₁(0)/√(1 − V²) as in 54.
- **63B (resting spin, static).** 59's resting polarised source packet (n = 97, as built by `exgc59.source_terms`). Magnetic moment m = ½ Σ x × s of the D2-E1 source against the PS source j = P + ¼Cᵀs, both spins.
- **Recompute.** R_dyn from the saved momentum histories by an independent least-squares fit; 63B moment by direct summation in a separate script. Shared: packet builders, operators, evolve.

## Predictions (written before any run)

Strength-to-R line from 52 (λ = 0, R = 1.1152) and 53 (λ = 4, R = 0.9975): R(λ) = 1.1152 − 0.02943 λ (54 confirmed linearity with λ = 2).

| | Prediction | Basis |
|---|---|---|
| P1 | in-place strength 4 Σs / Σ P at t = 0, pair: 3.85 ± 0.02 (equal to 54's uniform-N value) | flat weights: the operator no longer sees N |
| P2 | R_dyn(0.2) = 1.002 ± 0.008 | line at 3.85 (54's registered band) |
| P3 | \|R_dyn(0.2) − 1\| < 2% | Lorentz clause |
| P4 | F₁(0) within 1% of 53's 4.699 × 10⁻⁴ | the static pull is set by N, not the shift |
| P5 | 63B: \|m(D2-E1)\| / \|m(PS)\| < 0.05, both spins | in the continuum ½{Ĥ, α̂} = p̂ exactly ({α_i, α_j} = 2δ_ij, {β, α} = 0): the energy current has no spin-curl part, and a resting packet has no orbital current |

If P5 holds, the directed count supplies the orbital (moving-mass) shift but not the spin part, and 58's plaquette-circulation term ¼Cᵀs stays a separate counted thing. If P2 fails while P1 holds, the remaining gap is in the operator's shape, not its strength.

## Checks

- **C1.** Zero packet gives zero source (exactly).
- **C2.** D2-E1 source flips sign with V (V = ±0.2 at t = 0, static, relative difference < 1e-12).
- **C3.** With N ≡ 1 everywhere, D2-E1 equals 54's D2-E to 1e-14.

## Budget

54's D2-E at one V took 0.45 core-h (pair 1066 s, single 550 s). 63A: about 0.9 core-h. 63B and checks: under 0.05 core-h. Cap 2 core-h. Platform recorded with results.
