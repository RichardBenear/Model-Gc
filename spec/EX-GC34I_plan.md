# EX-GC34I — a lattice kinetic term derived from an exact re-slicing identity (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "Go ahead with deriving a kinetic term from an exact lattice
re-slicing identity". This plan and the `explore/LOG.md` entry are pushed before any 34I computation. Cap: 2
core-hours on the cloud; I expect under 0.2. Our construction, not LaViolette's Model G.

## Idea (derivation)

**Variables and action.** Per vertex: lapse n, shift β (3), and 7 Kuhn edge changes e. Edge strain rates are
σ_t = l⁰_t(ė_t − D̂_t·(β(v + d_t) − β(v))). The kinetic energy is T = ½ σ† 𝒦(q) σ for a 7×7 Hermitian symbol 𝒦(q).
The potential is the unchanged 34D lapse × Regge term, whose lapse–edge coupling is C_en(q).

**Two gauge symmetries.**
- Displacement symmetry, δe = D̂·Δu with δβ = u̇, is exact for any 𝒦, because σ is invariant (34E).
- Linearised re-slicing by a vertex scalar ξ is δn = ξ̇, δβ = μ g(q) ξ, δe = 0, where g is a lattice gradient
  and μ a normalisation.

**The identity.** In the continuum, re-slicing is a symmetry because the DeWitt contraction of ∂_i∂_jξ with ḣ is
ξ times the time derivative of the linearised scalar curvature. The lattice analogue follows by working through the
11×11 kernel: re-slicing is an exact null vector of M(q, ω) for every ω if and only if

  **(I)  𝒦(q) σ_ξ(q) = C_en(q) / l⁰**  (elementwise in the 7 edges), where σ_ξ = −μ l⁰ D̂·(e^{iq·d} − 1) g(q).

The shift rows then vanish automatically, because the Regge coupling is displacement-gauge invariant (34D).

**The derivation.** I look for local symbols 𝒦(q) that satisfy (I) exactly:
- entries are trigonometric polynomials with real-space offsets in a box of range R;
- real coefficients;
- Hermitian, i.e. K_{−o} = K_oᵀ.

(I) is linear in (𝒦's coefficients, μ). It is imposed as an exact polynomial identity on an FFT grid large enough
to contain every product offset. Among all solutions, I take the one closest to the 34E per-tetrahedron term T, in
the Frobenius norm of the real-space coefficients. For that comparison, T's own μ is set by least squares. Because the
target is T, the continuum limit of T is kept as nearly as (I) allows. If the solution space is empty, (I) has no
local solution at that range and g, which is a valid result.

## Candidates (fixed now; no others are scored)

| id | g(q) | range R |
|---|---|---|
| If1 | forward difference, g_a = (e^{iq_a} − 1)/h | 1 |
| Ic1 | central difference, g_a = i sin(q_a)/h | 1 |
| If2 | forward | 2 |
| Ic2 | central | 2 |

Range 2 is built and scored only for a g whose range-1 system has no exact solution.

## Criteria (thresholds fixed now)

- **I0: identity exact.** The constructed 𝒦 satisfies:
  - (I) to ≤ 1e-12 relative, at 200 random q;
  - Hermiticity to ≤ 1e-14;
  - M(q, ω) v_ξ = 0 to ≤ 1e-12 relative, for v_ξ = (−iω, μg, 0) at 50 random (q, ω).
  - Internal check: the 7×7 symbol extracted from 34E's T reproduces `exgc34g.kin_T` to ≤ 1e-13.
- **G1–G3, G5 and G7** as in `spec/EX-GC34G_plan.md`, with the same direction set, spacings, thresholds and
  argument-principle confirmation (`exgc34f.analyse_posthoc`), with one change. In G3 the expected number of exact
  re-slicing directions is now one in every direction, and the count of static Jordan chains must equal
  nullity_full11 − 3 as before.
  - G1: no growing root, Im ω > 1e-8 |ω|.
  - G2: at kh ≤ 0.05, exactly two confirmed real roots with 0 < ωh < 0.5, each within 1e-3 of the light speed.
    Every other non-zero root has ωh ≥ 0.5 and is reported with its kh scaling.
  - G5: static kernel equals 34D; gauge null; E4 gravitomagnetism ≤ 1e-2 at kh = 0.05.
  - G7: direct real-space sum of the kinetic energy on a periodic 8³ lattice matches the Fourier form to ≤ 1e-10.
- **G4: time domain.** Scored only for a candidate that passes I0, G1–G3 and G5. Generic and (1,1,0) directions at
  kh = 0.05, with the zero-shift-gauge descriptor system evolved as in 34G's G4 and the same thresholds.
  - The 34G G4 recipe assumed a regular edge kinetic block. If A_ee is singular for the candidate, G4 is reported as
    not applicable.

**A candidate passes** if it meets I0, G1–G5 and G7.

**My expectations, stated before running:**
- Range 1 may well have no exact solution. I'd guess range 2 does, for at least one g.
- If an exact solution exists, I expect G1 to pass, because the 34F growth came from the near-gauge re-slicing
  direction, which becomes exact gauge.
- G2 is the main risk. The 7th edge keeps kinetic energy here, so a third mode exists. Whether it sits at the cutoff
  (pass) or low (fail) I can't predict.

## Reported, not scored

- The solution-space dimension at each range.
- μ.
- The size of the correction ‖𝒦 − 𝒦_T‖ / ‖𝒦_T‖.
- The real-space support of the correction.
- The 34E wave-speed table.
- The third mode's frequency and kh scaling.

## Scope

Linear, flat background, continuous time. An exact linear identity does not guarantee a non-linear version. A
passing candidate would still need non-linear and coupled-matter tests.

## Files

`src/exgc34i.py`:
- `derive` → `results/EX-GC34I/derive.json` (the solutions and I0);
- `spec` → `candidates.json` (G1–G3, G5);
- `evolve` → `evolve.json` (G4);
- `check` → `check.json` (G7).

It reuses `src/exgc34e.py`, `exgc34f.py` and `exgc34g.py` unchanged.
