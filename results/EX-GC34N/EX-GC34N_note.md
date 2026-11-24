# EX-GC34N — a staggered lattice (Y1) makes both gauge symmetries exact and gravitomagnetism second order (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "I sign off on 34N").

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "close 34M and then redesign the shift couplings and vacuum
action".
- The plan `spec/EX-GC34N_plan.md` and the `explore/LOG.md` entry were pushed before any run (ef80f66).
- Code: `src/exgc34n.py`.
- Outputs: `build.json`, `spec.json`, `spec_mp.json`, `static.json`, `sync.json`, `matter.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**The redesign works.** One registered gate missed, narrowly, at the coarsest spacing.

**The new lattice.** Y1 places each field by the parity of its spatial indices:
- lapse, scalar matter and h_aa at sites;
- shift and displacement on links;
- shear h_ab on plaquettes.

Every gauge transformation is then a centred one-cell difference.

**The vacuum action** is linearised Einstein–Hilbert, with each derivative replaced by the matching staggered
difference.

**What it gives:**

| | J1 (34J–34L) | Y1 (34N) |
|---|---|---|
| re-slicing and displacement identities | both exact | both exact (5e-17) |
| locality | range 1 | range 1 (periodic to 8e-16) |
| growing roots | none | none (50-digit roots, ≤ 5e-13) |
| waves | 2 at c (7.5e-4 at kh 0.05) and a 7th-edge cutoff mode | exactly 2, ω = \|K\| (7e-5 at kh 0.05); no extra mode |
| static extra zero modes | none | none (nullity exactly 4 = static gauge) |
| static γ | 1 | 1 (2e-15) at kh = 0.05, 0.5, 2 |
| **gravitomagnetism E4 at kh = 0.05** | 1.25e-2, ratio 2.00 (first order) | **1.04e-4, ratio 4.00 (second order)** |
| synchronous time domain (G4″) | pass (34K) | pass |
| scalar packet, energy identity | exact (34M) | exact (≤ 9e-15) |
| **scalar packet, local momentum identity** | **first order, 14–20% at h = 0.25** | **0.2–0.4% at h = 0.25, about second order** |

**The registered miss.** E-M2 asks for the momentum-identity residual to fall by at least 3.5× per halving of h at
both steps. From h = 0.5 to 0.25 it does in every case (3.54–4.04). From h = 1 to 0.5 it does in only 3 of 6 cases;
the others are 2.87–3.44. At h = 1 the packet is only 3 cells wide (σ = 3, k₀h = 0.6), so that step is
pre-asymptotic. So **E-M2 FAIL as registered**, with second-order convergence on the finer step.

**Why the momentum fix worked.** In 34M the momentum density had to sit on links while J1's displacement acted at
vertices, which gave an O(1) first-order mismatch. In Y1 both live on links, and the mismatch drops by a factor of
50–70 at h = 0.25.

**What changed physically.** The spatial metric is no longer 7 Kuhn edge lengths but 6 staggered components:
normal components at sites, shear on plaquettes. The non-metric 7th-edge mode is gone. The lapse still multiplies
the linearised scalar curvature, so the 34C "pure clock" and single-covariant-energy structure that gave γ = 1 is
kept.

**Gravity and the massless scalar now share the same lattice dispersion**, ω = |K|, with K_a = 2 sin(q_a/2)/h. So
waves and matter see the same cone exactly, as 34A wanted. Along the axis it matches 34A's light speed exactly. Along
diagonals it is more isotropic than 34A's reference: 0.986 against 0.959 at kh = 1 on (1,1,1).

## Registered results

**A. Construction**

| | result |
|---|---|
| I0: max null residual over 20 random (q, ω) | 5.4e-17 |
| I0: Hermitian | yes |
| **I0** | **PASS** |
| G7: periodicity M(q + 2π e_a) = M(q) | 8e-16 |
| G7: stencil range | 1 |
| G7: largest entry outside range 2 | 1e-16 |
| **G7** | **PASS** |

**B. Spectrum** (synchronous gauge, 300 cases: 34F/34G directions at kh 0.0125–0.2, plus axis, (1,1,0) and (1,1,1)
up to the zone boundary)

| | result |
|---|---|
| G1, 50-digit roots: max \|Re λ\|/\|K\| | 4.8e-13 → **PASS** |
| G2: wave roots per case | exactly 4 (two ± pairs) |
| G2: zero cluster per case | 8 |
| G2: other roots | none |
| G2: max wave speed − 34A light speed at kh ≤ 0.05 | 6.9e-5 → **PASS** |
| G3: static nullity at every case | exactly 4 → **PASS** |

The zero cluster is the residual synchronous gauge: 4 of its 8 directions satisfy the constraints, as in 34J/34K.

**C. Static field and gravitomagnetism**

| | result |
|---|---|
| γ at kh = 0.05, 0.5, 2 (3 directions) | 1 to 2.2e-15 → **G5-static PASS** |
| E4 at kh = 0.0125 / 0.025 / 0.05 / 0.1 | 6.5e-6 / 2.6e-5 / **1.04e-4** / 4.17e-4 |
| E4 order ratio | **4.00** → **G5-E4 PASS** |
| Riemann change under each lattice gauge | ≤ 1.2e-16 (exactly invariant readout) |
| forced-equation residual (source exactly conserved on the lattice) | 1.5e-15 |

**D. Time domain (G4″)** — generic and (1,1,0), 40 periods, two steps:

| | result |
|---|---|
| constraints | ≤ 8e-13 |
| wave tidal amplitude change | ≤ 8e-13 |
| non-zero seeds | 1.9–2.7e-8 of the wave |
| zero-sector seeds | ≤ 5.4e-13 |
| steps agree | yes |
| **G4″** | **PASS** |

**E. Matter: the 34M packet on Y1**

| | result |
|---|---|
| E-M0 (Richardson) | 4e-12 / 6e-12 → **PASS** |
| E-M1, energy identity | ≤ 8.8e-15 → **PASS** |
| total momentum | conserved to ≤ 8e-14 |
| E-M2: r_P at h = 1 / 0.5 / 0.25, axis | 0.038–0.051 / 0.012–0.014 / 0.0030–0.0037 |
| E-M2: r_P at h = 1 / 0.5 / 0.25, generic | 0.026–0.028 / 0.0068–0.0082 / 0.0018–0.0022 |
| E-M2 ratios, coarse step (h = 1 → 0.5) | 2.87–3.78 |
| E-M2 ratios, fine step (h = 0.5 → 0.25) | 3.54–4.04 |
| **E-M2** | **FAIL as registered** (3 of 6 coarse-step ratios below 3.5) |

**Overall: not a pass as registered (E-M2).** Every other gate passes.

**My stated expectations:**
- Everything passes, E4 about 1e-4: right, 1.04e-4.
- E-M2 SECOND: the convergence is second order on the finer step, but the coarse step fell short of the registered
  3.5.
- My named risk, a first-order momentum mismatch in the off-axis shear coupling, did not occur.

## Recomputation and shared dependencies

- **Gauge invariance** is checked directly, with the null vectors applied to the assembled kernel. The kernel is built
  from the linearised Einstein tensor `G_up`. That is written independently of `exgc34e.gr_solution`, which builds
  its own G_of.
- **E4 compares two separate solves:** the lattice solve with K and `exgc34e.gr_solution` with k. Shared dependency:
  `exgc34e.riemann`.
- **Roots** were computed twice: float64 (`spec.json`) and 50-digit mpmath built from scratch (`spec_mp.json`). They
  agree on the wave roots. Only the float64 run shows the splitting artefact (see Deviations).
- **The matter sources** were checked against a directly coded Lagrangian (E-M0).

## What this means

- **The 34J–34M chain is resolved at the level of linear theory.** With every field on its index-parity position, both
  gauge symmetries are exact and local at once. The vacuum is second order. A dynamical matter field's sources
  conserve energy exactly and momentum to second order.
- **This is a different lattice from Regge/Kuhn.** The spatial metric is carried by staggered components, not edge
  lengths. That matters if Model Gc should keep "link weights = spatial metric" literally. Here the shear is on
  plaquettes and the normal components at sites.
- **Natural next steps (none approved):**
  1. Re-run E-M2 with a wider packet or finer spacings (h = 0.5, 0.25, 0.125) to settle the order as a registered
     result.
  2. The ChatGPT review's step 2 on Y1: couple the packet's own sources to Y1 and test whether energy changes gravity.
  3. Return to the knot on Y1.
  4. A white-paper update once that holds.

## Compute

About 1.3 CPU-minutes in total:
- build: 1 s;
- spec: 1.5 s;
- spec_mp: 46 s;
- static: 0.1 s;
- sync: 1.3 s;
- matter: 12 s.

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18, mpmath.

## Deviations

- **G1/G2 root computation.**
  - The first spec run used float64 eigenvalues. The 8-fold zero cluster of the synchronous system is exactly
    defective (Jordan blocks), so float64 rounding split it into spurious roots of size about (1e-16)^{1/4}, up to
    8.5e-5 |K|. As first computed that showed "G1 FAIL" and "other roots" (`spec.json`).
  - I rebuilt the system in 50-digit arithmetic (`spec_mp`). The split falls to 5e-13, consistent with
    (1e-50)^{1/4}, with exactly 4 wave and 8 zero roots everywhere.
  - I score G1 and G2 on the 50-digit roots, because the registered criterion concerns the roots, not a particular
    eigensolver. Both outputs are committed for your judgement.
- None otherwise. E-M2 is scored strictly as registered.
