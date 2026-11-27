# EX-GC48 — a Dirac-type lattice fermion under the 44 coupling rules

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC48_plan.md) (with amendment A1, made before any run) | about 0.03 core-h on the cloud container (Linux, Python 3.11, numpy). Our construction; not LaViolette's Model G, not Chain Theory.

## Conclusion

**One Wilson–Dirac fermion fixes both problems from 47 at once.** Its inertia equals its energy over c², up to a correction of order mass/cutoff. At rest it is dust: its spatial ("stress") charge goes to zero as the packet widens. When it moves it carries the relativistic stress, in proportion to v². Here c is the speed of the same field when it is massless, measured from a moving packet.

**But a composite of two such fermions bound by an on-site contact is not relativistic.** Its inertia is only 28–79% of its energy. Its uncancelled binding stress is 40–320 times its binding energy, against 0 in GR.

The reason is that a contact interaction on a lattice binds at the cutoff scale: the internal kinetic and potential energies are each 2–3, while m ≤ 0.5. The potential energy of an instantaneous site interaction carries no inertia and no stress.

**So the matter field is no longer the problem; the binding interaction is.** The interaction would have to be carried by a field that has its own momentum and stress. That leads into 49, a propagating firing-rate field.

**The 47 "inertia failure" stands under an independent signal speed.** The fastest group velocity of the 45/47 network is 0.174. With that, 2 t_pair E/v_max² = 0.0039 (B1) and 0.0042 (B2), against 1 for relativistic matter.

## Scorecard

| Prediction | Result | Score |
|---|---|---|
| P1 free R = 2 t_pair E_rest/c² = 1 + r m (1e-3), axis and diagonal | Agrees to < 1e-6 at all 10 (m, r) points, both directions | PASS |
| P1 massless c = 1 within 1% | Centroid speeds 0.9795 / 0.9948 / 0.9987 at k₀ = 0.4 / 0.2 / 0.1; extrapolated to k₀ → 0: 0.99999 | PASS (deviation D1) |
| P1-2 R₄₇ < 0.05 (47 failure survives) | 0.0039 (B1), 0.0042 (B2); v_max = 0.174 | PASS |
| P2 Q_N/E = 1 (1e-12) | 1 to rounding | PASS |
| P2 rest slope −2 ± 0.1 (σ 16→32) and within 5% of 3(1 + r m/2)/(2σ²m²) at σ = 32 | m = 0.5: −1.95, −1.1% (PASS). m = 0.2: −1.79, −5.4% (FAIL) | **1 of 2** |
| P2m moving: ΔQ_a/E = Δ⟨v²⟩(1 + r m/2)/(1 + r m)², within 5% | 1.0–3.6% at all six points | PASS |
| P2-W2 rest Q_a/E ≈ −3r/m | −14.4 (pred −15), −5.95 (pred −6) | PASS (W2 is not dust) |
| P3a Q_N = E₀ | Exact (homogeneous weights) | PASS |
| P3b δR < 0, within 30% of V/⟨E₁+E₂⟩ | δR < 0 at all 4 points; within 12%, 29%, 14%, but 34% at m = 0.2, E_b/2m = 0.02 | **3 of 4** |
| P3c η > 0.5 | 40–318 | PASS |
| P3d spectator formula within 5% of t_c | m = 0.5: 3.0%, 3.9% (PASS); m = 0.2: 7.4%, 8.4% (FAIL) | **2 of 4** |
| C1 L = 8 dense vs fast solver, 1e-9 | Energies 3e-14 (PASS); Q_a 1.2e-9 (FAIL, finite-difference noise floor at h = 1e-5) | **1 of 2** |
| C2 real-space ⟨H_link⟩ vs k-space Q_a, 1e-10 | 4e-16 (W1), 3e-15 (W2); total energy agrees too | PASS |

## Evidence

**48A, free inertia** ([A.json](A.json)). The positive band is E² = (t sin k)² + (m + t r W)². Its low-momentum curvature gives inertia m/(1 + r m), so R = 1 + r m exactly. Fermions much lighter than the cutoff 1/r are relativistic, with a correction of order m r. This is the reverse of the 45–47 design, where the rest mode had to be far above the bandwidth.

The signal speed was measured from massless packets. **Deviation D1:** the planned k₀ = 0.05, σ = 20 packet straddles k = 0, where the massless bands are degenerate and part of the packet moves left. I used k₀ = 0.4, 0.2, 0.1 with momentum width k₀/4 instead, and extrapolated in k₀².

**48A-2, the 47 bodies** (Codex's calibration). The largest Hellmann–Feynman group velocity over all seven bands of the 45/47 network is 0.1741 (B1) and 0.1740 (B2). The bodies' 2 t_pair E is 1.2e-4, so R₄₇ ≈ 0.004. With the network's own fastest signal, the pairs' inertia is about 250 times too large. Measured against the 37B/44 optical speed of 1, it is about 8,000 times.

**48B, one packet** ([B.json](B.json), [B_posthoc.json](B_posthoc.json)). Q_a/E at rest falls as σ⁻² toward 0, so a wide packet at rest is dust.

The m = 0.2 failure is a prediction error. I scored the asymptotic formula at σ m = 6.4, where the k²/m² corrections are still about 5%. Post hoc (not scored), at σ = 64 and 128 the slope is −1.98 and the formula holds to 0.4% (m = 0.2) and 0.07% (m = 0.5).

A moving packet's stress is v² times the lattice factor (1 + r m/2)/(1 + r m)². That factor is 0.78 at m = 0.2 and 0.56 at m = 0.5, and tends to 1 (the continuum's p·v) as m → 0.

Under the W2 reading, where the Wilson on-site part does not stretch, Q_a/E = −3r/m. The Wilson term must be treated as a whole lattice derivative.

**48C, the composite** ([C.json](C.json)). Values are at r = 1, W1, along the axis.

| m | E_b/2m | L | U | T | V | R_c | δR | V/⟨E₁+E₂⟩ | η | Q_a/E₀ |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.02 | 96 | 5.78 | 1.99 | −2.01 | 0.787 | −0.475 | −0.672 | 110 | 2.24 |
| 0.5 | 0.08 | 96 | 5.91 | 2.91 | −2.99 | 0.488 | −0.675 | −0.764 | 40 | 3.46 |
| 0.2 | 0.02 | 128 | 6.19 | 2.40 | −2.41 | 0.519 | −0.567 | −0.860 | 318 | 6.49 |
| 0.2 | 0.08 | 128 | 6.24 | 3.27 | −3.30 | 0.276 | −0.770 | −0.900 | 108 | 9.39 |

- The diagonal direction agrees with the axis to 0.3% at m = 0.5 and to 6% at the weakest m = 0.2 point.
- At 0.75 L the weakest point is not size-converged: δR moves from −0.519 to −0.567 and η from 282 to 318. The other three points move by less than 2%.
- A contact interaction needs U ≈ 6, about the cutoff, to bind at all in 3D. The bound state is then made of lattice-scale momenta, so the internal motion is ultra-relativistic.
- The composite's inertia comes mostly from the kinetic energy of its parts. The potential energy, about −2 to −3, adds none. Gravity sees the full E₀ through Q_N. So under the 44 rules the composite would fall 20–70% slower than a free fermion. This is an order-1 violation of the equivalence principle that depends on composition.

**Checks** ([checks.json](checks.json)). Two separate paths confirm the solver and the packet charges.
- C1 builds the L = 8 two-body Hamiltonian in the band basis from numerically diagonalised spinors. It reproduces the 5 × 5 projector solver's energies to 3e-14.
- C2 builds the real-space Wilson–Dirac Hamiltonian at L = 24. Its ⟨H_link⟩ for a positive-band packet equals the k-space Q_a to 3e-15, for both Wilson readings.
- Shared dependencies: numpy. C1 calls the fast solver for comparison and C2 calls the analytic Q_a formula.

## What this does and does not show

- It shows that relativistic single-particle matter can be put on the lattice with the 44 clock and space rules intact: inertia = energy, and dust at rest. The rest energy does not need to sit on a site.
- It does not give relativistic composites. Any instantaneous lattice interaction carries no momentum and no stress, so its potential energy has no inertia and its binding stress is not cancelled. A shallow, longer-range potential would shrink δR to about V/E, of order the binding fraction, but would not remove it.
- **Suggested next (not approved):** bind two fermions by exchange of a propagating lattice field (for example a massless or light scalar coupled to ψ̄ψ, with its own link and site terms under the 44 rules), and check whether the field's momentum and stress restore R_c = 1 and η = 0. This is a lattice von Laue test, a field version of idea 6.
