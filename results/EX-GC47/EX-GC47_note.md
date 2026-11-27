# EX-GC47 — two bound heavy pairs: they attract as predicted, but fall about 10⁴ times too slowly

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC47_plan.md), pushed before any run. PI request and go-ahead 2026-11-27 04:35 UTC ("the review of 46 is available. I would like to see the tests outlined there done"). Implements the three stages of the [Codex review of GC46](../../review/EX-GC46AB_codex_review.md). PI CLOSED 2026-11-27. Cloud container (1 core), Python 3.13, numpy 2.5. About 0.05 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**Two bound heavy pairs attract through the prescribed fields with exactly the expected coefficient, and the attraction moves them. But their acceleration is only about 1.2 × 10⁻⁴ of the Newtonian value, because a pair's inertia (1/(2 t_pair) ≈ 8,000) is not tied to its energy (≈ 1).** That was predicted before the run, and it is the main result. Two registered clauses fail (below).
- **Static attraction (47B).** At L = 32 and separations 3–12, the cross energy of two microscopic pairs equals the point-charge prediction −G(E₁E₂ + E₁Q_a2 + Q_a1E₂) in the same periodic box to 0.14%. The spatial cross terms add 1.68%, as Codex estimated (1.66%). It is reciprocal to 1e-13 and linear in G to 1e-13.
- **Release (47C).** A probe pair released in the source's field accelerates inward at −2 t_pair ∂V₁₂/∂r to within 0.8–2.6%. That is 1.16–1.27 × 10⁻⁴ of the acceleration a slow test body would have in the same lapse field. Released together, two pairs approach each other.
- **What it means.** Making the rest energy on-site (45/46) makes the body nearly stress-free, but it also removes the link between rest energy and inertia. Inertial mass ≠ gravitational mass, by a factor of about 8,000.

## 47A — benchmark pairs (L = 32, heavy band)

| | B1: (m₀, g) = (0.5, 0.025) | B2: (0.7, 0.025) |
|---|---:|---:|
| E_b/W_h | 0.15 | 0.15 |
| E (pair energy) | 0.99141 | 1.38951 |
| E_b | 1.121e-3 | 1.167e-3 |
| Q_a (t∂E/∂t) | 8.216e-3 | 1.453e-2 |
| Q_a/E | 0.83% | 1.05% |
| t_pair (fit, K → 0) | 5.936e-5 | 4.626e-5 |
| t_pair, 46 single step (four bands) | 5.881e-5 | 4.581e-5 |
| t spread over x, y, z, (1,1,1) | 2e-4 | 5e-5 |
| rms relative separation | 0.50 sites | 0.45 sites |
| probability both on one site | 78% | 81% |
| min breakup gap over sampled K | 1.121e-3 (at K = 0) | 1.167e-3 (at K = 0) |
| **m_i/E = 1/(2 t_pair E)** | **8,496** | **7,778** |

**Scoring P1: PASS.**
- t_pair agrees with 46's single-step value to 0.9% and 1.0% (gate 10%).
- Isotropy: spread 2e-4 and 5e-5 (gate 5%). This is guaranteed for the quadratic term by cubic symmetry, so it checks the code more than the physics.
- The pair is compact: rms separation 0.50 and 0.45 sites (gate 5).
- m_i/E = 8,496 and 7,778 (gate > 10³).

## 47B — static cross energy from the pairs' microscopic charge profiles

The pair's lapse and strain densities come from transition form factors F(q) = ⟨Ψ_q|o(q)|Ψ₀⟩ in the natural gauge (the band-projected resolvent of an on-site pair, which is gauge-invariant in the Bloch phases). Fields come from the 44 equations with the local vacuum counterterm, 4πG = 0.01, a periodic box, and the zero mode removed. The cross energy is the probe's first-order energy in the source's fields. The comparators are point charges (E, Q_a) in the same periodic Green function, and dust (E only).

| L | V₁₂/V_point − 1, 6 ≤ r ≤ 12 (all orderings, σ = 0 and 2) | V₁₂ < 0 for | first r with V₁₂ ≥ 0 (point comparator changes sign there too) |
|---|---:|---|---:|
| 16 | up to 9.4% (near the sign change) | r ≤ 6.9 | 6.93 |
| 24 | up to 9.0% (near the sign change); 0.28% for r ≤ 9 | r ≤ 10.4 | 10.39 |
| **32 (registered domain)** | **≤ 0.14%** | all r ≤ 12 | 13.86 |

- **Dust ratio (B1/B1, L = 32, x axis).** V₁₂/V_dust = 1.0168–1.0170 at r = 4–11. The prediction 1 + 2Q_a/E is 1.01671, using F_a(0).
- **Short range.** V₁₂/V_point at r = 2 is 1.0016 and at r = √3 it is 0.9993. The pair-size (form-factor) correction is below 0.2% everywhere tested, because the pairs are compact.
- **One-body control.** One heavy carrier as source and probe: V/V_point − 1 = 8e-6 at r = 6.
- **Box size (P4, descriptive).** V₁₂·r/G_eff at r = 6 is −0.08, −0.33 and −0.49 at L = 16, 24 and 32, against the bare −0.999. The absolute potential is dominated by the periodic zero-mode offset and has not converged at L = 32. The ratio to the point-charge comparator converges, and so does the force (47C uses gradients).

**Scoring P2.**
- **V₁₂ < 0: PASS on the registered domain** (L = 32, r = 3–12). At L = 16 and 24 V₁₂ becomes positive beyond r ≈ 6.9 and 10.4, exactly where the point-charge periodic potential does. That is the neutralising background of a periodic box, not the matter.
- **Within 3% of the point-charge prediction at r ≥ 6: PASS at L = 32** (≤ 0.14%). It fails at L = 16 and 24 next to the periodic sign change.
- **Dust coefficient 1 + 2Q_a/E within 3%: PASS** (≤ 0.14%).
- **Reciprocity to 1e-10: PASS** (≤ 1.5e-13). This is an algebraic identity of the 44 source equations; the test confirms the code respects it.
- **Linearity to 1e-9: PASS** (≤ 3.9e-13). V₁₂ = 0 exactly at G = 0.
- **F_N(0) = E to 1e-8: PASS** (≤ 1e-15).
- **F_a(0) = Q_a to 1e-6: FAIL.** F_a(0) = ⟨H_link⟩ = 8.282e-3, but t∂E/∂t = 8.216e-3, a 0.8% difference. In the band-projected model the Hellmann–Feynman theorem for strain is broken, because the contact interaction's band form factors depend on t. The effect on V₁₂ is about 0.014%. The four-band t-derivative from 46 is 8.219e-3, so the band projection error in Q_a is under 1%.

## 47C — release

**Reduction and its validation.** The probe's centre of mass evolves under H_eff = E(K̂) + V₁₂(R̂), with E(K) the exact pair band and V₁₂ from the K = 0 form factors. At L = 6 this is compared with the exact two-body problem: both carriers in the heavy band, the full 46,656-state pair space, and the same static source field.

| level | exact shift | H_eff shift | rel. diff. | full first-order (K-dependent F) rel. diff. |
|---|---:|---:|---:|---:|
| 1 | −1.5227e-3 | −1.5186e-3 | 0.27% | 1.36% |
| 2 | −2.5447e-4 | −2.5243e-4 | 0.80% | 2.17% |
| 3–5 | −2.3186e-4 | −2.2740e-4 | 1.92% | 2.09% |

**P3 first clause: PASS** (≤ 1.9%, gate 3%). The first run of this check missed degenerate copies of the six-fold level (Lanczos with k = 5); see Deviations. The K-independent H_eff does slightly better than the full first-order projection, so the residual is second order in the field, not the K approximation.

**Probe release** (L = 24, source B1 fixed at the origin, probe packet σ = 3 at rest, a source-off run subtracted; T chosen so the probe moves about 0.05 site).

| probe, start | a measured | −2 t_pair⟨∂V₁₂/∂x⟩ | rel. | a_N = −⟨∂Φ_S/∂x⟩ | **a/a_N** | dt halving | field change across the pair / E_b |
|---|---:|---:|---:|---:|---:|---:|---:|
| B1 at 6 | −1.788e-9 | −1.773e-9 | +0.8% | −1.497e-5 | **1.19e-4** | 6e-9 | 0.019 |
| B1 at 8 | −1.209e-9 | −1.233e-9 | −2.0% | −1.040e-5 | **1.16e-4** | 1e-7 | 0.009 |
| B2 at 6 | −1.892e-9 | −1.942e-9 | −2.6% | −1.497e-5 | **1.26e-4** | 1e-8 | 0.023 |
| B2 at 8 | −1.322e-9 | −1.350e-9 | −2.1% | −1.040e-5 | **1.27e-4** | 5e-8 | 0.011 |

**Both released** (relative coordinate, total momentum 0, start at 8). The separation falls by 0.044 (B1+B1) and 0.045 (B1+B2) sites in the run. The accelerations match −2(t₁ + t₂)⟨∂V/∂r⟩ to 3.7% and 0.02%. Equal and opposite momenta hold by construction in this coordinate.

**Scoring P3.**
- **H_eff vs exact within 3%: PASS** (1.9%).
- **Acceleration equals −2 t_pair ∂V₁₂/∂r within 5%: PASS** (0.8–2.6%).
- **a/a_N below 10⁻³: PASS** (1.16–1.27 × 10⁻⁴).
- **a(B1)/a(B2) differs from 1 by more than 20%: FAIL.** The ratio is 0.945 at 6 and 0.915 at 8, so free fall differs between the two bodies by 5.5–8.5%, not by 20%. The reason: a/a_N ≈ 2 t_pair E_P, and t_pair·E happens to be close for B1 (5.9e-5 × 0.99) and B2 (4.6e-5 × 1.39). The bodies still fall differently, but less than I predicted.
- **Both released: PASS.** The separation decreases in both cases, and dt halving changes the displacement by ≤ 2e-7 (gate 1e-3).

## Checks

- **C1:** F_N(0) = E PASS; F_a(0) = Q_a FAIL (0.8%, above).
- **C2 (cubic fill): PASS.** The wedge-filled F equals a direct evaluation at six random full-grid q to ≤ 9e-16.
- **C3 (exact vs H_eff): PASS** (above, after the eigensolver fix).
- **C4 (norm): PASS.** Norm is conserved to ≤ 4e-13. Energy conservation was not recorded (split-step with static H; not separately checked).
- **Separate recomputation (`src/exgc47_recompute.py`): PASS.** It uses an opposite bond convention, dense pair problems at all 216 K, its own gauge rule, real-space two-body amplitudes and real-space operators, with no symmetry fill. At L = 6 it reproduces F_N and F_a over all q to ≤ 1.4e-13, and V₁₂ at (3,0,0) and (2,1,0) to ≤ 2.4e-12. Shared: the model definition, U and the kernel formula.

## Deviations

- **C3 eigensolver.** The first run used eigsh(k = 5, ncv = 40) and returned only some copies of the six-fold degenerate COM level, so its no-field check was off by 6.2e-5 and its fifth level was misaligned (12% "error"). Rerun with k = 12, ncv = 96 (no-field check 2.9e-14); the first run is kept in C.json as `validation_first_run`. The dynamics were not affected.
- **Scoring domain for P2.** The plan put the main separations at L = 32 (3–12) and smaller boxes at shorter range. Clauses are scored on L = 32; the L = 16 and 24 failures near the periodic sign change are reported above.
- **Q_a convention.** 47B uses F_a(0) = ⟨H_link⟩ (what the fields see); 47A reports t∂E/∂t. They differ by 0.8% (C1).
- **Packets.** Release used σ = 3 only, two separations, L = 24. The F tables for 47C (`partial_F24.npz`, 112 kB, SHA-256 142d5b4de5e6d4c73869a7de4c621f1fa28714bd65c16c56cb9074cfaf2312fb) are regenerated by part B and not committed.

## What it means (inference)

- **The coupling works as built.** Microscopic bound bodies, with their own charge profiles, attract through the 44 fields with the predicted coefficient and reciprocity, and the attraction moves them. This was put in through the field equations, so it is a consistency result, not a derivation.
- **The equivalence of inertial and gravitational mass is missing.** A pair's gravitational charge is its energy (plus a 1–2% strain term). Its inertia is 1/(2 t_pair), set by how easily it hops, which is about 8,000 times larger. So it falls at about 10⁻⁴ of the Newtonian rate. In relativity m_i = E/c² comes from Lorentz invariance: a massive particle's dispersion is E = √(m²c⁴ + p²c²). Here the rest energy sits on a site and the mobility comes from weak hybridisation, so the two are unrelated. The route 45–46 used to make matter dust-like is exactly what breaks m_i = E.
- **What would fix both at once (suggestion, not run).** Matter with a relativistic lattice dispersion, such as a Dirac-type (Wilson) fermion whose mass term couples chiral components on a site and whose hopping obeys the same 1/a strain law. Then E(k) ≈ √(m² + c²k²) near k = 0. This gives m_i = E/c² automatically. It also gives Q_a = t∂E/∂t = c²k²/E, which is exactly zero at rest, so the body is dust-like with no small parameter. The natural next test is whether such a fermion's bound pairs keep m_i = E and have zero binding stress.

## Correction after the Codex review (added 2026-11-27; scoring above is unchanged)

The [Codex review of GC47](../../review/EX-GC47_codex_review.md) found two issues. Both are fixed in [EX-GC47R](../EX-GC47R/EX-GC47R_note.md).
- **Release packet seam.** The launch Gaussian was built from coordinates wrapped about the origin, so at L = 24, start 8, σ = 3, a site holding 61% of peak weight was given a displacement of −20 instead of +4. Rebuilt periodically about the packet's own centre and repeated over L, σ, start and fit window: acceleration within 4.4% of −2 t_pair∂V/∂r, a/a_N = 1.15–1.32 × 10⁻⁴, B1/B2 = 0.90–0.95. The conclusions above stand.
- **F_a(0) ≠ Q_a (C1).** The 0.8% gap was the projector-derivative term of the band-projected model. Adding it reproduces t∂E/∂t to 4e-10. V₁₂ changes by ≤ 0.018%.
- **Wording.** "Inertia ≈ 8,000 × energy" is in lattice units with hopping 1. A relativistic comparison needs the independently measured signal speed: the clean statement is 2 t_pair E ≠ c_signal².
