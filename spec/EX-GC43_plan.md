# EX-GC43 — what does the lapse couple to? One energy functional, positive-energy excitations, and the vacuum term

2026-11-26 | EXPLORATORY | Claude | Plan only. PI request in thread: "yes, write 43" (2026-11-26 18:19 UTC). PI run go-ahead 2026-11-26 18:23 UTC ("go ahead"), given after this plan was pushed and before any run. Follows the [Codex review of GC42](../review/EX-GC42_codex_review.md). Cap 2 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why

42 fed an *assigned* energy density into the lapse constraint. If the constraint comes from varying one functional with respect to N, the lapse couples instead to the conjugate density q_v = ∂H/∂N_v (Codex, "main gap").

Under 37B weights, H[N] = Σ N_vN_w H_link,vw + Σ N_v H_site,v. At N ≡ 1 this gives q_v = 2·(link energy at v) + (site energy at v). The total charge is therefore Q = ΔE + ΔE_link, not ΔE.

From the 42A data, the charge-to-energy ratio Q/ΔE₀ = d ln ΔE₀/d ln N is already non-universal for hubs: it runs 0.10–0.39 across the six hubs. That finding is post hoc, so here it becomes a scored prediction under a new method.

Two questions follow:
1. Under one consistent functional, is the gravitational charge of *free positive-energy excitations* of the unchanged network universal, Q/ΔE equal for all?
2. Can the vacuum's own term be removed by anything internal to the model, or only by a counterterm?

**Frame (inference, written before any run).**
- Link energy counts twice. That is the lattice version of Tolman's ρ + 3p: radiation has active mass 2E.
- In GR, a static bound body still has active mass E, because its stresses integrate to zero (von Laue).
- So a pinned hub, held by an external coupling, is not expected to satisfy Q = E. A free packet is also not expected to, because it is not stationary-bound.
- The open question is whether any state of this network has Q = ΔE.

## Model

- **Network.** As in 42: the 40A pair network with 37B weights. λ₀ = 0.1 (main) and 0.2 (repeat).
- **Functional.** F[N, Ψ] = ⟨Ψ|H[N]|Ψ⟩ + (1/8πG) Σ_links (n_v − n_w)² + C[N], with n = N − 1.
- **Constraint.** Stationarity in N gives Δn = 4πG (q_v + ∂C/∂N_v), with q_v = ⟨Ψ|∂H/∂N_v|Ψ⟩.
- **Prescriptions.** The counterterm C is chosen in one of three ways, all fixed now.
  - **V (raw variational).** C is the uniform cosmological term only, chosen to cancel the uniform vacuum tadpole σ_vac = ∂ε_vac/∂N at N = 1. This is put in by hand and the note will say so. At linear order the vacuum kernel is χ(q) = ∂q/∂n at wavevector q.
  - **L (local counterterm).** C = −Σ_v ε_vac(N_v), the local vacuum energy at uniform lapse N_v. The kernel becomes χ(q) − χ(0). This is 42's variant A, now derived from the functional.
  - **X (excitations only).** The source is q_v[Ψ] − q_v[vacuum of the same H[N]]. This subtraction is non-local. At linear order there is no vacuum kernel, and the far field is pure Poisson with charge Q.
- **Consequence.** At linear order, Q is the same under L and X. They differ only in G_eff. So universality is decided entirely by Q/ΔE.

## Parts

### 43A — vacuum response of the conjugate source

- σ_vac(λ₀) and χ(0) = ε_vac''(1), from uniform N = 1 ± 1e-3 (Bloch, L = 32).
- χ(q) from the 42B supercell: n = ε cos(qx), L = 32, m = 1..8, ε = ±1e-4. Measure q_v = Re Σ_{i∈v} (∂h/∂N_v · P)_ii, using the corrected q = π normalisation.
- Then the V and L far fields of a unit point charge, on 128³ as in 42C.

### 43B — box resonance (Codex item 3, descriptive)

- Under V, with the measured χ: the point-charge response n(r = 10) for periodic boxes Lg ∈ {96, 112, 128, 144, 160}.
- A Dirichlet comparison: the radial local-kernel equation (1 + 4πG c₂χ) ∇²n − 4πG χ(0) n = 4πG Q δ, with n(R) = 0 for R ∈ {30, 40, 50, 60}.

### 43C — hub charges under the conjugate source

- For the six 42 hubs, q_v from δP by Woodbury (as in 42's hub_source, with the ∂h/∂N_v weights).
- The integrated charge Q = Σ_v δq_v, compared against dΔE₀/dN at uniform N (central difference at N = 1 ± 1e-3).

### 43D — positive-energy excitations of the unchanged network

These are one added particle, or one particle–hole pair, in Gaussian packets of Bloch states, with no change to H.

- **Packets.**
  - Conduction band (lowest band above μ), centred at k = 0 and at the band minimum, with widths σ = 1, 2, 4 sites.
  - The same widths in the highest band.
  - A neutral particle–hole pair: a conduction packet plus a valence hole packet, same centre, σ = 2.
- **Charge.** The primary charge is the time-averaged (diagonal-ensemble) one, Q = Σ_k |c_k|² (2e_link(b,k) + e_site(b,k)). It is conserved under the free evolution, which answers Codex's "static or time window". The t = 0 frozen value is recorded as descriptive.
- **Energy.** ΔE = Σ_k |c_k|² e(b,k) for an added particle (canonical, no μ term, so that ΔE and Q are counted the same way). For a pair it is e_c − e_v, with Q = Q_c − Q_v. The ratio Q/ΔE is independent of normalisation.
- **Band survey.** Also the band-resolved ratio 1 + e_link/e over the full zone of every band, as a map. Here e_link = ∂e/∂t at t = 1, with t the hopping scale (Hellmann–Feynman).

## Predictions (scored)

- **P1.** χ(0) < 0 at both λ₀, and χ(0) ≤ 2 ε_link,vac.
  - Reason: ε_vac(N) = N² ε_link + N ε_site + (state relaxation). The relaxation is ≤ 0 at second order.
  - So the raw variational vacuum (V) also turns the constraint into a Helmholtz equation, and V gives no 1/r: the 1/r fit residual is > 10× that of L.
- **P2.** Under L, n is 1/r outside the core, with periodic-Green's-function fit residual < 1% at r = 4–20. G_eff/G = 1/(1 + 4πG c₂χ), with c₂χ from χ(q) − χ(0) ≈ c₂χ q̃². The sign of c₂χ is recorded, not predicted.
- **P3 (hubs).**
  - Q/ΔE₀ matches d ln ΔE₀/d ln N within 1%.
  - It is non-universal: the spread over the six hubs at fixed λ₀ exceeds 10%.
  - Q < 0, the same sign as ΔE₀.
- **P4 (free excitations: the main test).**
  - Every packet has Q > 0 and ΔE > 0: positive energy and positive (attractive) charge.
  - Q/ΔE > 1 for conduction-band packets, because their link energy is positive (mostly antibonding). The value is recorded; no upper bound is predicted, since e_site can be negative.
  - Q/ΔE is **not** universal: the spread over the packets exceeds 5% at fixed λ₀.
  - If the spread is below 1%, universality holds for free excitations. That would be a surprise, and the note would say so.
- **P5 (descriptive).** The band map of 1 + e_link/e: its range, and whether any band or region has ratio exactly 1, meaning energy that is all on-site.
- **P6 (box resonance, descriptive).** Under V:
  - n(r = 10) changes by more than 2× or changes sign across the Lg sweep.
  - The Dirichlet amplitude diverges near k_χ R = jπ, with k_χ² = 4πG|χ(0)|/(1 + 4πG c₂χ).

## Checks

- **C1 (Hellmann–Feynman).** At L = 8 (dense), q_v from P equals the finite difference of E with respect to a one-site N_v perturbation (ε = 1e-4), within 1e-6 relative, both for the vacuum and for one hub.
- **C2 (kernel limit).** The supercell χ(q → 0), after the q̃² trend is extrapolated out, equals the uniform χ(0) within 1%.
- **C3 (packet charge).** For one σ = 2 conduction packet at L = 16, the band formula equals the real-space dense computation of ⟨∂H/∂N_v⟩ summed over v, within 1e-10. The time average equals the mean over 200 evolution times within 1e-3.
- **C4 (hub charge).** The Woodbury Σ_v δq_v equals the uniform dΔE₀/dN within 1e-4 relative at λ₀ = 0.1. The tolerance is 1e-3 at λ₀ = 0.2, because of the cube truncation seen in 42's C1.

## What would change the picture

- **If P4 holds,** the consistent coupling makes mass depend on how energy is stored. That is the same N-scaling defect behind 42A's equivalence-principle failure, so it would close the "pair network plus lapse" line unless some class of states has a universal ratio.
- **If any class does (P5),** for example states whose energy is all link energy (ratio 2, "built like light"), the next study builds composites from that class only, and asks whether a bound composite recovers Q = E through a lattice von Laue balance.

## Outputs

`src/exgc43.py`, `results/EX-GC43/`. A separate recomputation path, `src/exgc43_recompute.py`, recomputes Q/ΔE from the saved band data, independently of the driver's packet code.
