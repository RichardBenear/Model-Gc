# EX-GC44 — separating the clock from space: lapse homogeneity, a spatial link factor, and whether a bound composite is stress-free

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC44_plan.md), pushed before the PI's run go-ahead ("go ahead and execute the plan for doing those", 2026-11-27 00:00 UTC) and before any run. It follows the [Codex review of GC43](../../review/EX-GC43_codex_review.md). PI CLOSED 2026-11-26. Cloud container (1 core), Python 3.13, numpy 2.5, scipy. Under 0.1 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**Separating time from space fixes the equivalence-principle failure, and shows exactly what is left over.**
- **Weight the whole Hamiltonian by the lapse (44B), and mass equals energy, identically.**
  - Every hub and every packet has Q_N = ΔE, and the redshift test passes to 1e-15.
  - But light's speed then scales only as N, so it bends at the Newtonian value, half of GR's.
- **Give space its own link factor (44C), and the stress Q_a = ΔE_link becomes a second, separate charge.**
  - It is the lattice version of 3p. 43's charge is recovered exactly as Q_N + Q_a.
  - For single excitations, Q_a ≈ 0.8–1.25 ΔE, so light from them bends 1.41–1.63× more than light from "dust". This matter behaves like radiation.
- **Binding does not cancel the stress (44D).**
  - A mobile exciton exists, and binding lowers its link energy. But the excess stress vanishes only if the attraction *grows* as links stretch: w* = −0.13 to −0.16 at λ₀ = 0.1, and −0.49 to −0.52 at λ₀ = 0.2.
  - At λ₀ = 0.1, w* drifts 16% with binding strength. At λ₀ = 0.2 it drifts only 5%, which is below the registered 10% test, so that clause fails.
  - The constituents themselves carry stress at rest. The unbound pair at the band edges has link energy 2.0, against energy 1.4 or 0.8.
- **The end-to-end chain holds (44A).** For localized sources, the far-field amplitude follows the conjugate charge to 2e-6 (by Gauss's law).
  - At equal energy (λ₀ = 0.2, ΔE = 0.90006), a neutral pair sources 55.3% more lapse than a charged particle.
  - The registered Green's-function fit missed its 1e-3 gate for the anisotropic neutral sources.

## 44B — lapse-homogeneous weights (control)

Links are weighted N̄ = (N_v + N_w)/2, and sites are weighted N. At uniform N, H = N·H(1).

- **Charge.** Q_N/ΔE = 1 for all 7 packet types and all 6 hubs at both λ₀. The largest deviation is 3e-12 (packets) and 1e-9 (hubs, finite difference).
- **Redshift.** ΔE₀(1.1)/1.1 ÷ ΔE₀(0.9)/0.9 − 1 is at most 4e-15 for every hub.
- **Optics.** The top-band group velocity scales as N^p, with p = 1.0000000000000 at both λ₀.

**P2: PASS.** This rule is 42's variant-A assigned density: q_v = site + ½ Σ links = ρ_v. So 42's universal 1/r was this homogeneous model all along (inference, by inspection).

## 44C — two fields: Q_N (energy) and Q_a (strain stress)

Hop ∝ N̄/a. Q_a = −a ∂H/∂a at a = 1, by uniform finite difference.

| State (λ₀ = 0.1 / 0.2) | ΔE | Q_N/ΔE | Q_a/ΔE | Light bending vs dust, 1 + Q_a/(2Q_N) |
|---|---|---:|---|---|
| Flat conduction, any σ | 0.900 / 0.800 | 1 | 1.111 / 1.250 | 1.556 / 1.625 |
| Top band σ = 1 | 0.987 / 1.037 | 1 | 0.988 / 0.815 | 1.494 / 1.407 |
| Top band σ = 2 | 0.925 / 0.879 | 1 | 1.070 / 1.042 | 1.535 / 1.521 |
| Top band σ = 4 | 0.907 / 0.822 | 1 | 1.100 / 1.180 | 1.550 / 1.590 |
| Particle–hole σ = 2 | 1.425 / 0.879 | 1 | 1.397 / 2.179 | 1.698 / 2.090 |
| Hubs (6) | negative | 1 | −0.68 to −0.90 / −0.59 to −0.79 | — |

**P3: PASS.**
- Q_N = ΔE (within 1e-8), and Q_a reproduces 43's ΔE_link.
- The single-packet bending ratio is 1.41–1.63, inside the predicted 1.4–1.7.
- No state has Q_a = 0. Hubs have Q_a of the opposite sign to their energy: the defect is held under tension by its imposed coupling.
- C0 (unit-source field solves): the Ψ and Φ coefficients equal the analytic values to 2e-16.

## 44D — the exciton

The exciton is a particle in the top band and a hole in the top valence band, with an on-site density attraction −U. It is solved exactly in the band-projected two-body problem, using its separable (rank-36) form, at L = 32. Binding energies are checked at L = 24 and L = 12.

| λ₀ | E_b / gap | U | ⟨V⟩ | δ⟨H_link⟩ (bound − band edges) | w* = −δlink/⟨V⟩ | E₀(K) − E₀(0), K = (1,0,0)·2π/32 |
|---|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.05 / 0.15 / 0.30 | 2.640 / 3.309 / 4.253 | −0.521 / −0.718 / −0.964 | −0.081 / −0.103 / −0.129 | **−0.156 / −0.144 / −0.133** | −8e-4 / −9e-4 / −1.2e-3 |
| 0.2 | 0.05 / 0.15 / 0.30 | 5.857 / 6.352 / 7.039 | −0.892 / −1.071 / −1.266 | −0.464 / −0.543 / −0.625 | **−0.520 / −0.507 / −0.493** | −3e-4 / −3e-4 / −4e-4 |

For w ∈ {0, 0.5, 1, 2} the excess stress δQ_a = δlink + w⟨V⟩ is negative in every case: −0.08 to −3.2. That is net tension.

**P4:**
- **A bound, mobile exciton exists: PASS.** It has nonzero dispersion. Its minimum is not at K = 0: E₀ falls slightly as K grows.
- **No lattice von Laue identity (w* changes by > 10%): PASS at λ₀ = 0.1 (16% change), FAIL at λ₀ = 0.2 (5.2%).** The registered "surprise" threshold of 2% is not reached either: w* is not constant.
- **Physically,** every w* is negative. Stress-free binding would need an interaction that strengthens as links lengthen (U ∝ a^{+0.13 to +0.52}), which is confinement-like. No ordinary (w ≥ 0) scaling cancels the stress.

**P5 (descriptive).** The unbound pair at the band edges has link energy 2.00 at both λ₀, against energy 1.40 and 0.80. So Q_a/ΔE = 1.43 and 2.50 before any binding. Binding removes only 4–8% (λ₀ = 0.1) or 23–31% (λ₀ = 0.2) of that.

## 44A — end-to-end

**Deviation.** The registered recipes were single-band packets: top band × valence, and flat band × valence. They are not spatially localized, so far-field fits of them are meaningless.
- At λ₀ = 0.1 (σ = 2), only 67% of the top-band packet's charge lies within r = 6, and 85% for the flat band.
- The cause is that the top band touches the flat bands at k = 0, and the band gauges have obstructions.
- I replaced them, before scoring, with packets built by projecting a fixed orbital (one link-end mode "e0", or the site-symmetric "uni") onto the gapped conduction or valence group. These are > 97.7% inside r = 6.
- No equal-energy pair of *neutral* sources with different link fractions is reachable with these orbitals. Their energy windows do not overlap, so that clause is untested. The equal-energy charged-vs-neutral case (Codex's second case) was built at λ₀ = 0.2.

| λ₀ | Source | ΔE | Q | Q/ΔE | Registered fit: amplitude ÷ (G_eff Q) | Gauss law, R = 10 / 16 / 24 |
|---|---|---:|---:|---:|---:|---|
| 0.1 | neutral e0 × uni | 1.4328 | 3.4203 | 2.387 | 1.0145 | 0.999991 / 0.999992 / 0.999993 |
| 0.1 | neutral e0 × e0 | 1.8168 | 3.8105 | 2.097 | 1.0269 | 0.999991 / 0.999992 / 0.999993 |
| 0.1 | charged uni σ = 2 | 0.9398 | 1.9253 | 2.049 | 0.9991 | 0.999993 (all R) |
| 0.2 | neutral e0 × uni | 0.90006 | 2.79223 | 3.102 | 1.0212 | 0.999975 / 0.999974 / 0.999976 |
| 0.2 | neutral e0 × e0 | 1.6527 | 3.5967 | 2.176 | 1.0315 | 0.999975 / 0.999974 / 0.999976 |
| 0.2 | charged uni σ = 2 | 0.9135 | 1.8020 | 1.973 | 0.9989 | 0.999976 / 0.999976 / 0.999977 |
| 0.2 | charged uni, σ tuned to ΔE = 0.90006 | 0.90006 | 1.79739 | 1.997 | 0.9985 | 0.999976 / 0.999976 / 0.999977 |

**P1:**
- **Amplitude ratio = charge ratio within 1e-3, with the registered periodic-Green's-function fit: FAIL** for the neutral sources (1.5–3.2% high).
  - Their profiles are anisotropic, because the e0 orbital is a +x link end. The axis-ratio check at λ₀ = 0.2 shows x-axis excesses of 12% at r = 4, falling to 4% at r = 20, and y/z deficits of 4% and 0.1%.
  - The exact lattice Gauss-law estimator, added after this was seen (post hoc), gives amplitude ÷ (G_eff Q) constant to 2e-6 across all sources and radii. The 1e-5 offset from 1 is the c₂χ estimate.
- **Equal-energy neutral pair difference > 5%: not tested** (deviation above).
- **Charged vs neutral at equal energy reproduces 43's 56% within 1%: PASS.** Q ratio 1.5535 here vs 1.5571 in 43, with different states. By Gauss's law the amplitudes follow that ratio.

## Checks

- **C0: PASS** (2e-16).
- **C1: PASS.** The real-space profile sums equal the band formula to ≤ 1e-15.
- **C2 (mapping a = 1/N̄ gives 43's charge):** packets PASS (≤ 7e-12). Hubs reach 1e-9 to 1.2e-6, against the 1e-8 gate, so they FAIL as registered. The likely cause is the O(h²) truncation of the h = 1e-3 hub finite differences.
- **C3: PASS, method changed.** At U = 0 the continuum minimum equals the band-edge sum exactly. At L = 10, a dense diagonalization of the full projected matrix equals the separable solver to 6e-15. The planned real-space Lanczos was replaced by this dense solve.
- **C4: PASS.** ⟨V⟩ = E₀ − ⟨D⟩ equals the Hellmann–Feynman U ∂E₀/∂U to ≤ 4e-10.
- **Separate recomputation** (`src/exgc44_recompute.py`): a dense L = 12 projected matrix, with eigenvalue finite differences only. It reproduces w* within 0.04% (λ₀ = 0.1) and 1% (λ₀ = 0.2) of the L = 32 values, which is the size of the finite-size shift.

**Other deviations.**
- U was found by inverting the separable condition exactly, not by bisection.
- 44D ran at L = 32 with an L = 24 cross-check, not L = 16 dense.
- 44A uses static t = 0 group-projected sources. Their Q is not exactly time-conserved, because they mix bands, and no time window was computed.

## What it means (inference)

- **The split works.** Separating the clock (N, on everything) from space (a, on links) gives mass = energy for the lapse, universally and with the equivalence principle intact. That is the GR structure. Earlier, 37B forced the two together.
- **What remains is the matter.**
  - In this network every excitation, even at rest, keeps most of its energy in hopping. So its spatial charge, the stress Q_a, is comparable to its energy.
  - The bending of light past such matter, and the gravitational effect of the matter itself, both carry this radiation-like extra.
  - In GR, massive matter has negligible rest stress, and binding cancels the rest (von Laue).
  - Here neither holds. Binding reduces the stress only partly, and cancelling it would take a confinement-like interaction.
- **Next.** Matter whose rest energy is *on-site*, with Q_a = 0 at rest, is the missing ingredient. One example is a heavy on-site mode, a "mass term" that does not hop, with light-like hopping only for its motion.
  - The question is whether a composite of such constituents has zero excess stress, given ordinary interactions.
  - That is a well-posed next test, and it does not need a new substrate.

## Correction after review (2026-11-27, Codex review of GC44)

The registered P4 scoring above stands. The interpretation in the headline and in "What it means" was wrong on one point.
- w* cancels only the *binding* change of the stress. The gravitational question is whether the *whole* composite has Q_a = L_bound + w⟨V⟩ = 0.
- From the saved D.json, that requires w = +3.68 / +2.64 / +1.94 (λ₀ = 0.1) and +1.72 / +1.36 / +1.09 (λ₀ = 0.2) at binding fractions 5 / 15 / 30%.
- These are positive and strongly state-dependent. So "cancellation needs an interaction that strengthens as links stretch" applies to the binding increment only.
- What survives: no single strain law cancels the total charge across states. See [EX-GC45](../../spec/EX-GC45_plan.md).

## Files

[Code](../../src/exgc44.py) · [Recompute](../../src/exgc44_recompute.py) · `A.json` (44A), `BC.json` (44B, 44C), `C0.json`, `D.json` (44D), `recompute.json`
