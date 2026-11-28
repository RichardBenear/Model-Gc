# EX-GC64 — does the per-local-tick count transfer, is it one interaction, and does it coexist with the spin term?

2026-11-28 | EXPLORATORY | Claude | Plan only; not run. PI asked for it ("yes", 18:25 UTC, to "Should I write the 64 plan?"), responding to the Codex 63 review. Go-ahead to run pending. Requested cap: 2.5 core-h (expected about 1.9). Our construction; not LaViolette's Model G and not Chain Theory.

## Why

EX-GC63 counted 54's directed signal count in local ticks (rule D2-E1: source 4κ·½{Ĥ₁, v̂₁}, flat hop weights) and the co-moving pair came out at R_dyn(0.2) = 1.0028. The Codex 63 review (`review/EX-GC63_codex_review.md`) asked for three things:
1. **Transfer.** Freeze the rule and predict other boosts, a second width and other well depths, with a smaller time step. The depth change directly tests the clock-well explanation of 54's deficit.
2. **Reciprocity.** Matter feels the shift through P while the source is the counted current. EX-GC55 showed that this mixed pairing does not conserve energy (drift 4.2 × 10⁻³ of the interaction energy, against 4 × 10⁻⁵ for P/P and for reciprocal J_E/J_E).
3. **Hybrid spin.** Keep 58's plaquette term ¼Cᵀs and replace only the orbital current by the count. From 63's moments this hybrid has (0.022832 + 0.249385 − 0.024568)/0.249385 = **0.9930** of the PS moment. Score 59's precession against that number with nothing adjusted.

What this does not test: where the 4 comes from (it stays a supplied factor; the count is an operator, not a tally of events), or the dimension dependence.

## Rules

- **P** (53/55): source and coupling are the link momentum P.
- **D2-E1** (63, frozen): source 4κ × flat-weight energy current; coupling P. Source refreshed once per step and held across the four RK4 stages, as in 54/63.
- **RJE1** (new, reciprocal): 55's reciprocal energy-current interaction H_β = −Σ_a ½{b_a, Ô_a} with Ô built from flat-weight Ĥ₁, v̂₁; the same operator sources the shift. Source evaluated at every stage (55's evolve). At N ≡ 1 it is exactly 55's RJE.
- **Hybrid** (spin, static): j = PS − P + s(D2-E1) = s(D2-E1) + ¼Cᵀs.

## Parts

All pair runs use 55's smaller box: n = 49, separation 14, T = 8, dt = 0.1, spin up, κ₀ fixed as in 55, σ₀ = S.SIG = 5.

- **64A transfer.** D2-E1 and P at the settings below. F₁ is fitted as in 54/55. R = F₁(V)/F₁(0)/√(1 − V²).

  | Setting | κ | σ | V |
  |---|---|---|---|
  | base | κ₀ | 5 | 0.2 |
  | shallow | κ₀/2 | 5 | 0.2 |
  | deep | 2κ₀ | 5 | 0.2 |
  | narrow | κ₀ | 4 | 0.2 |
  | slow | κ₀ | 5 | 0.1 |
  | reversed | κ₀ | 5 | −0.2 |
  | half step | κ₀ | 5 | 0.2, dt = 0.05 (D2-E1 only, with its own V = 0 at dt = 0.05) |

  **F₁(0):** at V = 0 both sources vanish at t = 0, so each (κ, σ) setting's F₁(0) comes from one P run. This is checked once (C2).
- **64B in-place strength (static, t = 0).** s = 4 Σ source_x / Σ P_x in the pair for D2-E1, and for 54's N-weighted D2-E (reported) at κ₀/2, κ₀, 2κ₀ and σ = 4.
- **64C reciprocity.** RJE1 in the base pair (V = 0.2; F₁(0) from P as above). Energy bookkeeping is not rerun: at N ≡ 1, RJE1 = 55's RJE and D2-E1/P = 55's MIX, so 55's numbers apply exactly (checked, C3).
- **64D hybrid spin.** 59's static pipeline (n = 97, zero box, 48³ probe, PS probe coupling) with the hybrid source; axis 12, axis 16, plane 16, both source spins, against the PS source rerun in the same script.
- **Recompute.** F₁ by an independent least-squares design matrix and by a t ≥ 2 fit. The hybrid moment by direct summation. Shared: packet builders, operators, evolve routines, probe Hamiltonian.

## Predictions (written before any run)

ΔR ≡ R(rule) − R(P) at the same setting. Using the 52/53 slope dR/dλ = −0.02943 and an effective strength λ_eff:

ΔR = 0.02943 × (4 − λ_eff).

λ_eff is not fitted. It is fixed by the t = 0 in-place strength s: λ_eff = s for D2-E1, and (s/4)² × 4 for RJE1 (the same factor in source and coupling). The preferred-frame part of R scales as V², so ΔR(0.1) = ΔR(0.2)/4.

| | Prediction | Basis |
|---|---|---|
| P1 | D2-E1 in-place strength 3.85 ± 0.02 at κ₀/2, κ₀, 2κ₀ (independent of depth). At σ = 4 it is reported, not gated (expected slightly lower: a narrower packet samples more lattice dispersion) | flat weights do not see N |
| P2 | D2-E1 ΔR = 0.02943 (4 − s) ± 0.003 at base, shallow, deep, narrow, reversed (base value +0.0044) | 63 line; band from Codex's fit-window spread |
| P3 | slow: ΔR(0.1) = ΔR(0.2)/4 ± 0.0015 | V² scaling |
| P4 | half step: R(D2-E1) changes by < 0.001 | source held per step is a time-step effect |
| P5 | RJE1 ΔR = 0.02943 (4 − s²/4) ± 0.004 (≈ +0.0086); \|R − 1\| < 2% | reciprocal: source and coupling both carry s/4 |
| P6 | hybrid / PS precession rate = 0.993 ± 0.003 at axis 12, axis 16 and plane 16, both spins | 63 moments; the difference is a compact orbital current, so higher multipoles are negligible at d ≥ 12 |

Reading:
- If P1–P3 hold, the per-local-tick count transfers across depth, width and boost: 54's deficit was the clock-well weighting.
- If P5 holds, the reciprocal per-local-tick interaction is consistent but carries λ_eff ≈ 3.71. One interaction then needs 4/(s/4)² ≈ 4.31 rather than 4, or a different "count × velocity" operator. That is the direct cost of reciprocity.
- If P6 holds, the counted orbital source coexists with the spin term at the measured normalisation.

## Checks

- **C1 (63's C2 revisited).** The velocity-flip sum with compensated summation (`math.fsum`), and the full source arrays of a mirror-reflected packet (x → −x, spinor transformed by the parity matrix the code's α_x convention needs) compared link by link. This decides whether 63's 1.5 × 10⁻¹¹ is rounding. Reported, not gated.
- **C2.** At the base setting, F₁(0) from D2-E1 equals F₁(0) from P to 10⁻³ relative.
- **C3.** At N ≡ 1, RJE1's H_β and source equal 55's RJE to 10⁻¹⁴, and D2-E1's source equals 55's MIX source.
- **C4.** The hybrid source's z-moment over PS's equals 0.9930 to 10⁻⁴ (identity from 63's numbers).

## Budget

The costs come from 55 (P pair + single at n = 49, T = 8: about 200 s per V; RJE about 720 s per V) and 63 (D2-E1 scaled to n = 49, T = 8: about 350 s per V).
- 64A: about 2900 s.
- dt/2: about 1100 s.
- 64C: about 720 s.
- 64D: about 400 s.
- Checks: under 100 s.

Total about 5300–6800 s (1.5–1.9 core-h). Requested cap 2.5 core-h. If the first few runs project past the cap, I stop and ask before the half-step runs. Platform recorded with results.
