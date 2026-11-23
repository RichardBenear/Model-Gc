# EX-GC34J — the re-slicing-exact kinetic term with continuum normalisation, first-order cure and synchronous-gauge time domain (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "you can close 34I. Go ahead with the suggested next step."
This plan and the `explore/LOG.md` entry are pushed before any 34J computation. Cap: 2 core-hours on the cloud; I
expect under 0.3. Our construction, not LaViolette's Model G.

## Background and a pre-registration probe

**From 34I.** The identity K σ_ξ = C_en/l⁰, with a forward-difference re-slicing gradient, has exact range-1
solutions. With μ = −1 (post hoc in 34I), the solution closest to the per-tetrahedron term T:
- has no growing root and exactly two isotropic waves at c;
- puts the 7th-edge mode at the cutoff;
- but has a first-order E4 error (1.25e-2 at kh = 0.05; T is second order, 4.2e-4).

**Probe (done before writing this plan, on 34I's post hoc kernel; not a 34J result).** Moving the moment source or the
shift readout by half a cell does not cure the first-order error. The best placement is the unshifted one,
1.25e-2. So the error sits in K itself, not in where things are placed.

## Candidates (fixed now; no others are scored)

All use the forward gradient g_a = (e^{iq_a} − 1)/h, **μ = −1 fixed in advance** (the continuum normalisation), and the
34I derivation code (`src/exgc34i.py`, `derive_one`), with real, Hermitian, range-R coefficients.

| id | R | selection within the exact-solution space |
|---|---|---|
| J1 | 1 | closest to T in the Frobenius norm of the real-space coefficients (34I post hoc If1m1, now registered) |
| J2 | 1 | first match T's first-order Taylor term Σ_o o ⊗ K_o = Σ_o o ⊗ T_o (least squares if inexact), then closest to T |
| J3 | 2 | as J2, at range 2 |

The rationale for J2 and J3: T is second-order accurate in E4. Matching its odd (first-order) part, where the identity
allows it, should remove the first-order error. Whether the identity allows it is part of the question.

## Criteria (thresholds fixed now)

- **I0: identity exact.** As in 34I: identity ≤ 1e-12, Hermiticity ≤ 1e-14, re-slicing null vector ≤ 1e-12.
  Reported: the first-order match residual for J2 and J3.
- **G1, G2, G3 and G7** exactly as registered in 34G and 34I, through the 34G machinery. That includes the same
  direction and spacing set, argument-principle confirmation, and G2 within 1e-3 of the light speed at kh ≤ 0.05.
- **G5.** Static kernel equals 34D, gauge null, and E4 ≤ 1e-2 at kh = 0.05 in all 9 cases.
  Added (scored): **E4 order.** The ratio E4(0.05)/E4(0.025) ≥ 3, meaning second-order convergence (an ideal ratio
  of 4). It is reported at kh = 0.0125, 0.025, 0.05 and 0.1.
- **G4′: synchronous-gauge time domain**, replacing 34G's descriptor recipe, which is ill-posed once re-slicing is
  exact (34I).
  - Gauge n = 0, β = 0. The system is the regular 7-edge ODE A_ee ë = C_ee e.
  - Generic and (1,1,0) directions at kh = 0.05. Implicit midpoint for 40 periods of the slower wave, at steps
    dt = 0.02/ω_max and dt/2.
  - Runs:
    - (a) the wave alone;
    - (b) the wave plus every eigenmode outside the zero cluster (|ω| ≤ 0.01 k), each seeded at 1e-8;
    - (c) the wave plus the zero-cluster (residual-gauge) invariant subspace, seeded at 1e-8.
  - The observable is the tidal field E_ij = ½ ḧ_ij, from the 34D metric map (h00 = 0 in this gauge).
  - Pass if, for both steps:
    - the lapse and momentum constraint residuals are ≤ 1e-9 relative;
    - the wave's tidal amplitude changes by ≤ 1e-6;
    - max_t |E_b − E_a| ≤ 1e-6 max_t |E_a|;
    - max_t |E_c − E_a| ≤ 1e-6 max_t |E_a|.
- **A candidate passes** if it meets I0, G1–G3, G5 (with E4 order), G4′ and G7.

**My expectations, stated before running:**
- J1 passes everything except E4 (first order, 1.25e-2) and E4 order.
- J2 or J3: I'm not sure the identity allows T's first-order term. My guess is that range 2 does and range 1 does
  not, and that matching it gives second order.
- I expect G4′ to pass where G1–G3 pass, because the synchronous system is regular. Run (c) tests that the residual
  gauge carries no tidal field.

## Reported, not scored

- Solution-space dimensions.
- ‖K − T‖/‖T‖.
- The third mode's frequency.
- E4 at four spacings.
- Real-space support.

## Scope

Linear, flat background, continuous time. A passing candidate would still need non-linear and coupled-matter tests,
and a closed-form understanding of the kinetic term.

## Files

`src/exgc34j.py`:
- `derive` → `results/EX-GC34J/derive.json` + `kernels.npz`;
- `spec` → `candidates.json`;
- `e4` → `e4.json`;
- `sync` → `sync.json` (G4′);
- `check` → `check.json`.

It reuses `src/exgc34e.py`, `exgc34f.py`, `exgc34g.py` and `exgc34i.py` unchanged.
