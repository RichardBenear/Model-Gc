# EX-GC34I — an exact lattice re-slicing identity gives a stable kinetic term; registered verdict fails on normalisation and first-order E4 (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "you can close 34I").

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "Go ahead with deriving a kinetic term from an exact lattice
re-slicing identity".
- The plan `spec/EX-GC34I_plan.md` and the `explore/LOG.md` entry were pushed before any run (a744f26).
- Code: `src/exgc34i.py`.
- Registered outputs: `derive` → `derive.json` + `kernels.npz`, `spec` → `candidates.json`, `check` → `check.json`.
- Post hoc outputs: `results/EX-GC34I/posthoc/`.

This is our construction, not LaViolette's Model G.

## Conclusion

**No candidate passes as registered.** The derivation itself works, and it gives the first lattice kinetic term in
this series with no growing mode at all.

**The derivation.** Requiring re-slicing to be an exact null vector of the 11×11 kernel for every frequency gives one
linear identity, K(q) σ_ξ(q) = C_en(q)/l⁰. Here C_en(q) is the lapse–edge coupling in the 34D lapse × Regge term,
σ_ξ(q) is the edge strain rate produced by a re-slicing, and l⁰ is the unperturbed edge length.
- With a forward-difference lattice gradient (g = (e^{iq} − 1)/h), the identity has exact solutions.
  - The coefficients reach only nearest-neighbour cubes (range 1).
  - Hermiticity holds.
  - The solution space is 123-dimensional.
- With a central-difference gradient, there is no local solution at range 1 or 2: the identity is missed by 29% and
  14% relative.

**The registered candidate If1** (forward gradient, range 1, the solution closest to the 34E per-tetrahedron term T):

| check | result |
|---|---|
| I0 identity exact | PASS (identity 6e-15, null vector 5e-16, Hermiticity 4e-16) |
| G1 no growing root | **PASS** (0 growing roots over the whole 34F/34G case set, 0 unconfirmed) |
| G2 exactly two low modes at c | **FAIL** (exactly two low modes in every case, isotropic, but both 1.1% slow: 0.9887 vs 0.9999) |
| G3 chains = exact symmetries | **PASS** |
| G5 benchmarks | **FAIL** (static and gauge exact; E4 gravitomagnetism 2.3e-2 > 1e-2) |
| G7 real-space recomputation | PASS (5.6e-16) |
| G4 time domain | not earned |

**Why G2 failed.** The registered rule fixed the re-slicing normalisation μ from a small-q least-squares fit to T, at
μ_T = −0.9776. But T does not satisfy the identity even at leading order: its best-fit residual is 20–35% at
q ~ 1e-3, direction-dependent. So that fit is not a continuum value. The wave speed comes out as √|μ|, which is
0.9887. This is a flaw in the plan's normalisation rule, not in the derivation.

**Expectations stated before running:**
- I guessed range 1 would have no solution: wrong for the forward gradient, right for the central one.
- I expected G1 to pass: held.
- I expected G2 to be the main risk because of a third mode: the third mode sits at the cutoff (ωh ≈ 6 for
  kh ≤ 0.1), as hoped. G2 failed only on the normalisation.

## Post hoc (labelled; after the registered run)

**1. μ = −1, the continuum normalisation** (candidate If1m1, same derivation, closest to T at μ = −1).
- G1, G3 and G7 pass.
- G2 passes: two modes everywhere, max |ω/k − c| = 7.5e-4 at kh ≤ 0.05.
- G5 still fails narrowly: E4 = 1.25e-2. The E4 error is first order in kh:

  | kh | E4 error |
  |---:|---:|
  | 0.0125 | 3.1e-3 |
  | 0.025 | 6.3e-3 |
  | 0.05 | 1.25e-2 |
  | 0.1 | 2.5e-2 |

  T was second order (4.2e-4 at kh = 0.05). The one-sided gradient makes the shift response first-order accurate.
  Putting β half a cell along each axis in the readout does not help (1.5e-2).

**2. The 34G time-domain recipe does not work here, for a structural reason.** With exact re-slicing, the zero-shift
gauge leaves a residual gauge freedom: a re-slicing combined with a compensating displacement. The (n, edge) pencil
is then singular at every ω, with a smallest singular value of about 1e-18 at arbitrary ω. A descriptor
time-stepper on a singular pencil is ill-posed, and round-off blew up to inf in both runs. The spurious descriptor
eigenvalue (ω/k = −2.61 + 0.17i) is lapse-dominated and appears at any ω. It is not a root found by the
argument-principle analysis.

**3. Synchronous gauge (n = 0, β = 0)** is available only because re-slicing is now exact. It gives a regular
7-edge system A_ee ë = C_ee e. Over the 34G case set plus the zone boundaries:
- Every case has four zero roots (residual gauge), two waves at c and one cutoff mode at ω ≈ 120 k (ωh ≈ 6).
- The only "growing" flags are split zero roots with |ω| ≤ 2e-3 k, which is numerical.
- The curvature of the zero sector equals the displacement-gauge control:
  - ratio about 1 in 183 cases;
  - up to 69 in two near-axis cases, where both values are below 5e-7;
  - at most 6e-7 where the control is exactly zero.
- Time domain at kh = 0.05, generic and (1,1,0), 40 periods, implicit midpoint:
  - constraints kept (lapse ≤ 2.6e-10, momentum ≤ 3e-12);
  - wave amplitude change 2e-8 to 7e-7;
  - the seeded-mode measure grew by a factor 2–7 and depended on the step. I think this is round-off exciting the
    residual-gauge linear drift and polluting the projection, but I have not shown it. So the time domain is
    clean for the wave and the constraints, but not yet clean on the seed measure.

## What this means

- **Physics.** Building in exact re-slicing removes the growing mode for structural reasons. The near-gauge direction
  that grew in 34E–34H is now exact gauge. The 7th (non-metric) edge combination keeps kinetic energy and sits at
  the cutoff instead of producing a slow third mode or a Jordan chain. That is the structure 34F–34H were missing.
- **Open defects:**
  - (a) The normalisation has to be fixed from the continuum (μ = −1), not from T. That is a rule change for a new
    registration.
  - (b) With the forward gradient, the shift response is only first-order accurate. A centred alternative has no
    local solution at range ≤ 2, so a staggered-shift formulation or a larger range would be the next step.
  - (c) A time-domain run in synchronous gauge needs a cleaner seed diagnostic, for example a tidal-field
    observable rather than a projection.
  - (d) It is linear only, and the kinetic term is a fit in coefficient space, closest to T, not a closed formula.
- **A sensible next study (not approved)** would register μ = −1 up front, compare forward with staggered β (or range
  3 central), and score the time domain in synchronous gauge with a tidal observable.

## Compute

About 3.5 CPU-minutes, well under the 2 core-hour cap:
- derive: 45 s;
- spec: 20 s;
- post hoc: about 2 minutes, mostly the blown-up descriptor runs.

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18.

## Deviations

- **First derive run.** It chose μ jointly with K, by minimising the distance to T (μ = −2.58). That did not
  match the plan's rule, so I corrected the code to the plan's rule (μ from T's least squares) and re-ran before
  scoring. The joint value is reported as `mu_joint_reported`.
- **Choice of q range for T's least squares.** I used small q (Gaussian, scale 0.02 per component), because T's identity residual is
  O(1) and direction-dependent. The plan did not specify the range.
- **Post hoc work.** Everything in the post hoc section was done after seeing the registered results.
