# EX-GC47R — corrections to EX-GC47 after the Codex review

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC47R_plan.md), pushed before any run. PI request and go-ahead 2026-11-27 13:32 UTC ("You probably should make the corrections it suggests"). Follows the [Codex review of GC47](../../review/EX-GC47_codex_review.md). PI closure pending. Cloud container (1 core), Python 3.13, numpy 2.5. About 1.2 core-h (most of it the L = 48 form-factor table). 47's registered scoring is unchanged; this is a corrected record. Our construction; not LaViolette's Model G, not Chain Theory.

**Both corrections leave 47's conclusions standing.** With the release packet built periodically about its own centre, and with Hellmann–Feynman-consistent form factors, the probe still accelerates at −2 t_pair ∂V/∂r (within 4.4% over 32 runs) and at 1.15–1.32 × 10⁻⁴ of the Newtonian rate. The microscopic force equals the point-charge force to 0.04% at every box size, and approaches 1/r² as the box grows. One registered clause fails: the energy-drift gate (1e-10) was too strict for the split-step integrator (worst 2.3e-7), though the trajectories are converged in dt to 1e-7.

## 1. Form factors from the derivative of the full projected operator

In the band-projected model the projector depends on the field. The local density the field couples to is therefore ⟨Ψ_q|P o(q) P|Ψ₀⟩ plus projector-derivative terms ⟨Ψ_q|V Q R o(q)|Ψ₀⟩ + ⟨Ψ_q|o(q) R Q V|Ψ₀⟩, per carrier, with single-carrier denominators.

| | B1, L = 16 / 32 | B2, L = 16 / 32 |
|---|---|---|
| 47's F_a(0) = ⟨H_link⟩ | 8.2824e-3 | 1.46746e-2 |
| correction, all other bands (variant i) | −6.685e-5 | −1.4223e-4 |
| **corrected F_a(0)** | **8.21555e-3** | **1.453234e-2** |
| t∂E/∂t of the projected pair | 8.21555e-3 | 1.453234e-2 |
| relative difference | 3.2e-10 / 3.6e-10 | 2.6e-10 / 2.6e-10 |
| correction, upper bands only (variant ii, Pauli-allowed) | −6.492e-5 | −1.4234e-4 |
| lapse correction at q = 0 | 4e-22 | (zero) |

- **Q1: PASS** (gate 1e-7). The 0.8% gap in 47 was the missing projector-derivative term, as Codex said. Restricting the intermediate states to the three other bands above μ changes the corrected value by 0.02% (B1) and 0.007% (B2).
- **Q2: PASS.** Using the corrected tables changes V₁₂ by at most 0.018% (L = 24) and 0.017% (L = 32) at r = 3–12 (gate 0.1%). Codex estimated 0.014%.

## 2. Release with the packet periodic about its own centre

The packet is now built from displacements wrapped about its own centre, and the centroid is tracked unwrapped about that centre. Previously, at L = 24, start 8, σ = 3, the site opposite the centre carried 61% of the peak weight with a wrong displacement. Now it carries at most 3e-4 (L = 24, σ = 3) and 7e-7 (L = 32, σ = 3). The sweep covers L = 24 and 32, σ = 2 and 3, starts 6 and 8, fit windows T and T/2, and dt and dt/2: 32 probe runs and 32 both-released runs.

| quantity | range over all runs | 47 |
|---|---|---|
| a_meas / (−2 t_pair⟨∂V/∂x⟩) − 1 | −4.4% … +3.4% | −2.6% … +0.8% |
| a/a_N, probe B1 | 1.15–1.20 × 10⁻⁴ | 1.16–1.19 × 10⁻⁴ |
| a/a_N, probe B2 | 1.24–1.32 × 10⁻⁴ | 1.26–1.27 × 10⁻⁴ |
| a(B1)/a(B2), same (L, start, σ, window) | 0.902–0.946 | 0.915–0.945 |
| dt halving, change in displacement | ≤ 1.0e-7 | ≤ 1.3e-7 |
| energy drift / |⟨V₁₂⟩| | 6e-12 … 2.3e-7 | not recorded |
| field-off centroid drift (subtracted) | ≤ 1.1e-3 site (L = 24, σ = 3); ≤ 4e-6 (L = 32) | |

**Scoring.**
- **Q3: acceleration within 5% of the prediction: PASS** (all 64 runs).
- **Q3: a/a_N within 10% of 47: PASS** (largest change +4%).
- **Q3: energy drift below 1e-10: FAIL.** The worst case is 2.3e-7, at σ = 3 with the full window. The Strang split-step conserves a shadow energy, not ⟨H_eff⟩ exactly. The gate was set too tight for this integrator. Since dt halving moves the trajectory by ≤ 1e-7, the accelerations are not affected.
- **Q4: B1/B2 ratio in 0.90–0.97: PASS** (0.902–0.946).

The residual spread in the acceleration ratio (±4%) comes from the fit window and packet width, not from dt. Shorter windows and σ = 3 agree best with the prediction (≤ 1.8%), because the packet then samples a smaller range of the band's non-parabolicity.

## 3. Force against distance and box size (B1 source and probe)

F(r) = −ΔV₁₂/Δr along x. The comparators are the point-charge force in the same periodic box (F_point) and the infinite-box force from −G_eff(E² + 2EQ_a)/r, with the same central difference (F_∞).

| r | F/F_point − 1 (L = 24 / 32 / 48) | F/F_∞, L = 24 | L = 32 | L = 48 |
|---|---|---:|---:|---:|
| 2 | 3.3e-3 (all L) | 1.105 | 1.106 | 1.107 |
| 3 | 2.5e-3 | 1.131 | 1.136 | 1.138 |
| 4 | 3.8e-4 | 1.059 | 1.071 | 1.077 |
| 5 | 2.4e-4 | 1.003 | 1.028 | 1.040 |
| 6 | 9e-5 | 0.951 | 0.997 | 1.019 |
| 7 | 5e-5 | 0.889 | 0.969 | 1.005 |
| 8 | 3e-5 | 0.806 | 0.936 | 0.993 |
| 10 | — | 0.531 | 0.843 | 0.966 |
| 12 | — | — | 0.685 | 0.928 |

- **Q5: PASS.** The microscopic force equals the periodic point-charge force to ≤ 3.8e-4 at r = 4–L/3, at every L (gate 0.5%). The pairs' internal size matters only at r ≤ 3 (0.3%).
- **Inverse square (descriptive).** At fixed r, F/F_∞ rises toward 1 as L grows. That is the periodic images' pull being removed. At L = 48 it is within 2% of 1/r² for r = 6–9. At r ≤ 4 the lattice Laplacian's Green function is 6–14% stronger than 1/r², independent of L: lattice discreteness, not the matter.

## What changes in 47's conclusions

Nothing substantive. The attraction is pure monopole far field with the predicted coefficient. The release acceleration tracks the force over the pair's inertia. The inertia-to-energy ratio of about 8,000 (acceleration about 1.2 × 10⁻⁴ of Newtonian) is unchanged. Codex notes that comparing an absolute mass with an energy needs the signal speed in the same units. The relevant test is 2 t_pair E_rest = c_signal², with c_signal measured independently in the light sector. That is the recommended next experiment, and it needs its own plan.

## Files

[Code](../../src/exgc47r.py) · `F.json` (Q1, Q2) · `R.json` (release sweep) · `force.json` · corrected F tables `partial_Fcorr24.npz` and `partial_Fcorr32.npz` are regenerated by part F and not committed.
