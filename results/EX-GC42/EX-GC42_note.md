# EX-GC42 — a pair network coupled to an energy-sourced lapse: does a hub source a universal 1/r?

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC42_plan.md), committed before any run. PI closure pending. Cloud container (1 core), Python 3.13, numpy 2.5. About 0.45 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**Yes, but only after the vacuum's own energy is subtracted by hand, and the "mass" comes out negative.**
- **Variant A (local vacuum energy subtracted): universal 1/r. PASS.**
  - Every hub (1 or 7 sites, Δλ = 0.1, 0.2, 0.4, both λ₀, both G) sources n = −G_eff ΔE₀ / r outside the core, with fit residual at most 6e-5.
  - G_eff/G is the same for all hubs to 3e-5, and equals 1/(1 + 4πG c₂) with c₂ ≈ 0.49 to 5e-6. The vacuum screens G slightly (0.5% at 4πG = 0.01).
- **Variant B (all energy gravitates): no Newtonian law.**
  - The vacuum's energy responds to N with K(0) = dρ_vac/dN ≈ −6.0. That turns Poisson into a Helmholtz equation, so the far field oscillates and grows into box modes. A 1/r fit is off by a factor of 2.6–19.
  - The registered k_J fit fails, because the periodic box resonates at k_J (see P2).
- **The hub's mass is negative.** A hub lowers the network's energy (ΔE₀ < 0), so it sources a lapse *above* 1 near it. In this network a hub is an anti-gravitating source relative to the vacuum.
- **The equivalence principle fails under 37B weights, as predicted.** The hub's proper energy ΔE₀/N changes by 11–17% between N = 0.9 and 1.1.
- **Honest framing.** At linear order, universality of the coefficient is guaranteed once the vacuum is gapped and its local response is subtracted. The source integrates to ΔE₀, and the kernel is the same for every hub. The run confirms that, and that the O(G²) corrections are small (0.07%, C2). It does not show that the subtraction itself is natural.

## 42A — uniform N (equivalence principle)

ρ_vac(N) scales as N^1.98 (λ₀ = 0.1) and N^1.91 (λ₀ = 0.2), because hopping (∝ N²) dominates.

| λ₀ | Δλ | ΔE₀ at N = 0.9 / 1.0 / 1.1 | Change of ΔE₀/N, 0.9 → 1.1 | Exponent p in ΔE₀/N ∝ N^p |
|---|---|---|---:|---:|
| 0.1 | 0.1 | −0.1056 / −0.1068 / −0.1077 | −16.5% | −0.90 |
| 0.1 | 0.2 | −0.2655 / −0.2709 / −0.2751 | −15.2% | −0.82 |
| 0.1 | 0.4 | −0.6976 / −0.7221 / −0.7428 | −12.9% | −0.69 |
| 0.2 | 0.1 | −0.1622 / −0.1662 / −0.1694 | −14.6% | −0.78 |
| 0.2 | 0.2 | −0.3652 / −0.3771 / −0.3868 | −13.3% | −0.71 |
| 0.2 | 0.4 | −0.8558 / −0.8936 / −0.9262 | −11.5% | −0.61 |

The gap changes too (1.08 / 1.40 / 1.76 at λ₀ = 0.1), and so do the in-gap levels: at λ₀ = 0.1, Δλ = 0.2, they move from 0.65 to 0.82 to 1.02.

**P1 (EP fails, > 1%): PASS (the expected failure).**
- ΔE₀ barely changes with N (it scales as about N^0.1–0.4), so the hub's proper energy ΔE₀/N is not a constant.
- A lattice-scale hub is not "built like light": its hopping-to-on-site ratio changes with N.

## 42B — vacuum response kernel K(q)

K is defined by δS_q = K(q) n_q, for the raw source (variant B). It was computed from an L = 32 supercell along x, for m = 1..16, with ε = ±1e-4.

| λ₀ | K(0) (uniform) | K at m = 1 | K(0) extrapolated from m = 1, 2 | c₂ = (K − K(0))/q̃² over m = 1 … 15 | K(π) |
|---|---:|---:|---:|---|---:|
| 0.1 | −6.00069 | −5.98167 | −6.00068 | 0.4949 → 0.4926 | −4.030 |
| 0.2 | −6.01020 | −5.99184 | −6.01019 | 0.4777 → 0.4715 | −4.125 |

- **P2, first clause (K(0) < 0 at both λ₀): PASS.** K(0) ≈ 2e_link + e_site, as Hellmann–Feynman says.
- **P4 (descriptive).** K(q) − K(0) is analytic: c₂ q̃², with c₂ constant to 0.5% (λ₀ = 0.1) and 1.3% (λ₀ = 0.2) across the whole zone.
  - G is renormalised by 1/(1 + 4πG c₂), with c₂ > 0, so the vacuum screens G slightly.
  - **Inference.** Most of c₂ is bookkeeping. Each site carries half of each link's energy, and the link weight is N_vN_w = 1 + n_v + n_w + …. That alone gives K − K(0) = −(e_link/2) q̃², which is about 0.5. The state's own polarisation is the small remainder.
- **Normalisation fix (recorded).** At m = 16 (q = π), cos² = 1 at every site, so the projection needs 1/L, not 2/L. The kernel run wrote 2/L there (−8.06, −8.25). The far-field driver halves it, and the code is fixed for future runs. Only the core region (q̃² > 4) depends on it.

## 42C — hub far field

S_hub is the hub's flat-space energy-density change from Woodbury, over a 13³ cube around the hub at L = 16. n is computed on a 128³ periodic grid, with K(q) interpolated in q̃² (held constant at the q = π value beyond the axis range).

**Variant A.** The fit is n = 4π (G_eff M) g on the axes at r = 4–20, with g the periodic lattice Green's function.

| λ₀ | 4πG | G_eff M / (G ΔE₀), range over 6 hubs | Spread | 1/(1 + 4πG c₂) | Largest fit residual |
|---|---:|---|---:|---:|---:|
| 0.1 | 0.01 | 0.995080 – 0.995086 | 6e-6 | 0.995079 | 2.4e-5 |
| 0.1 | 0.05 | 0.975868 – 0.975892 | 2.5e-5 | 0.975870 | 5.6e-5 |
| 0.2 | 0.01 | 0.995245 – 0.995259 | 1.3e-5 | 0.995255 | 2.7e-5 |
| 0.2 | 0.05 | 0.976679 – 0.976710 | 3.2e-5 | 0.976715 | 5.6e-5 |

**P3 (variant A: universal 1/r): PASS.** The residual is below 1% (it is below 6e-5), and the coefficient is the same for every hub within 1% (it is within 3.2e-5). Since ΔE₀ < 0, n > 0 near every hub.

**Variant B.**

| λ₀ | 4πG | √(4πG|K(0)|) | Fitted k (r·n = a cos(kr + φ), r = 4–40) | k ÷ k_J | 1/r fit residual |
|---|---:|---:|---:|---:|---:|
| 0.1 | 0.01 | 0.245 | 0.010 (lower end of the scan) | 0.04 | 19 |
| 0.1 | 0.05 | 0.548 | 0.499 | 0.91 | 2.6 |
| 0.2 | 0.01 | 0.245 | 0.010 (lower end of the scan) | 0.04 | 19 |
| 0.2 | 0.05 | 0.548 | 0.367 | 0.67 | 4.7 |

The same values hold for every hub, because they differ only by the overall factor ΔE₀.

**P2, second clause (oscillatory, k_J within 10% at both G): FAIL as registered.**
- The far field is not 1/r: the residual is 2.6–19 against the variant-A value of 6e-5.
- But it is not a clean cos(k_J r)/r either. With K(0) < 0, the constraint becomes ∇²n + k_J² n = 4πG S, a Helmholtz equation. In a periodic box, the box modes with q close to k_J resonate.
  - At 4πG = 0.01, the m = 5 mode of the 128 box has q = 0.2454, within 0.2% of k_J = 0.2450. That mode dominates: r·n first crosses zero near r ≈ 13 and then grows to 5× its near value.
  - At 4πG = 0.05, the fitted k is within 9% (λ₀ = 0.1) and 33% (λ₀ = 0.2).
- **Inference.** The deeper point is that a Helmholtz constraint has no unique decaying solution. A constraint with the raw vacuum as a source does not fix the lapse at all far from matter: the answer depends on the boundary.

## Checks

- **C1 (ΣS = ΔE₀ to 1e-8 at L = 16): PASS at λ₀ = 0.1, marginal FAIL at λ₀ = 0.2.**
  - At λ₀ = 0.1 the difference is ≤ 6.4e-11.
  - At λ₀ = 0.2 it is 1.1e-8 to 2.4e-7 (largest for the 7-site hub at Δλ = 0.2). The likely cause is the 13³ cube cutting off the longer λ₀ = 0.2 tail.
  - The far-field fit uses ΔE₀ from the Woodbury determinant, so the effect on the ratios is at most 1e-7.
- **C2 (dense self-consistent vs linear order, L = 8, 4πG = 0.01, within 5% at r = 2–4): PASS.**
  - The two agree to 0.065–0.066% (variant A) and 0.067–0.089% (variant B), after 28–29 mixing iterations (last change 7e-12).
  - The dense one-site kernel gives K(0) = Σ_r K(r), matching the Bloch derivative.
- **C3 (K(q → 0) → K(0) within 2%; diagonal vs axis within 5%): PASS.**
  - m = 1 is 0.32% from K(0) before extrapolation, and 6e-7 after.
  - Along (1,1,0), K agrees with the axis value at equal q̃² to 5e-7 (m = 1) and 2e-6 (m = 2).

- **Separate recomputation** (`src/exgc42_recompute.py`, from the stored axis profiles; it shares far_field and ΔE₀ but not the Green's-function fit). It fits n = −a/r + b + c r² on r = 6–30, a continuum 1/r plus the periodic box's constant and jellium terms.
  - The hub-to-hub spread is again below 2e-5 at every (λ₀, G), so universality is confirmed.
  - The absolute ratio comes out 1.15% higher than with the lattice Green's function (1.0065 instead of 0.9951 at 4πG = 0.01). That is expected: the continuum 1/r form leaves out the lattice Green's function's short-range corrections. The G_eff values in the table are from the lattice fit.
- The source file stores S to 12 significant digits, to keep the repository small. Recomputing the far field from it changes the ratios by 4e-13.

**Deviations.**
- The far field used a 128³ grid, not 64³.
- The kernel axis was extended to m = 1..16, to cover the whole axis range.
- The q = π normalisation was fixed in the analysis, as above.
- The variant-B fit used r = 4–40 along +x.

## What it means (inference)

- **A source law exists, but only with a tuned vacuum.** Coupling the pair network to a constraint-fixed lapse gives exactly what 40C lacked: a hub of any size and strength sources one universal 1/r, with mass equal to its energy and one slightly screened G.
  - This works only because the vacuum's local response ρ_vac(N_v) is removed by hand.
  - Let all energy gravitate (variant B), and the vacuum's own N-dependence acts as a negative mass term, which destroys Newton's law at the length 1/k_J = 1/√(4πG·6) ≈ 0.41/√(4πG) (4 sites at 4πG = 0.01). This is the cosmological-constant problem in miniature, and the network does not solve it.
- **The sign problem.** A hub that shares more entanglement lowers the energy, so it gravitates as negative mass. Positive gravitating mass needs something that costs energy above the vacuum: an excitation across the gap, or a hub with Δλ < 0. This is a pointed result for the "mass comes from shared entanglement" idea: in this network, more sharing means less energy.
- **The equivalence principle is still open.** Under 37B weights a lattice-scale hub's proper energy depends on N (exponent −0.6 to −0.9), so heavy and light hubs would not fall alike. 39B's "built like light" requirement is not met by on-site defects.
- **Next questions, in order of value:**
  1. Is there a natural reason the local vacuum term drops out? For example, a source defined as an energy *difference* that every event shares, or a lapse that counts excitations rather than energy.
  2. Does a positive-energy excitation (one particle across the gap) source the same G_eff?
  3. Can a hub be built whose ΔE₀/N is N-independent?

## Files

[Code](../../src/exgc42.py) · `ep.json` (42A), `kernel.json` (42B; the m = 16 entries carry the old 2/L normalisation), `source.json` (hub sources, C1), `far.json` (42C), `c2.json` (C2), `recompute.json` · [recompute](../../src/exgc42_recompute.py)
