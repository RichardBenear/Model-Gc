# EX-GC49 — firing-rate changes carried by the same coincidences as light: the second-order metric (PPN β)

2026-11-27 | EXPLORATORY | Claude | Plan only, pushed before any run. PI go-ahead 2026-11-27 15:52 UTC ("Let's do 7 then 1 then 4"; idea 1 of the GC40-to-now recap). Cap 2 core-h (expected under 0.3). Our construction; not LaViolette's Model G, not Chain Theory.

## Why

37B showed that event counting (clocks are one-site events at rate N, signals need two-site coincidences at rate N_xN_y) gives the weak-field metric g₀₀ = −N², g_ij = N⁻²δ_ij with γ = 1, but only for N given from outside. 37C showed that N must obey a wave equation with the light speed, or a moving mass has a preferred frame and drag. The open idea (37C, "what it means") is that a change in a site's firing rate is passed to its neighbours by the same two-site coincidences that carry light. Then gravity and light would share one cone automatically, and the field equation for N would come from the rule rather than be put in.

The first thing such a rule decides, beyond linear order, is the nonlinear completion of N(r) around a mass, i.e. the PPN parameter β in g₀₀ = −(1 − 2U + 2βU²). GR has β = 1. Perihelion advance is proportional to (2 + 2γ − β)/3, and solar-system tests bound |β − 1| below about 10⁻⁴. With γ = 1 already fixed by 37B, β is a sharp mechanism-level test.

## The rule

The 37B wave rule in each site's own time τ (dτ = N dt) is ∂_τ² u_x = Σ_y N_x N_y (u_y − u_x)/h². Idea 1 applies the same rule to the firing-rate field itself, with matter as a source s:

∂_τ² f(N)_x = Σ_y N_x N_y [f(N_y) − f(N_x)]/h² − s_x.

Here f(N) is the quantity a coincidence passes on. The literal reading of idea 1 is **f = N** (the firing rate itself). For a slowly varying N the right side is ∇·(N² f′(N) ∇N), so in a static vacuum F(N) = ∫₁ᴺ n² f′(n) dn is harmonic: F = −Q/r with Q = S/4π for a total source S (Gauss's law).

Four transported variables are registered, with their continuum predictions (N² = 1 − 2U + 2βU², U = Q/r):

| f(N) | Harmonic F | Exterior N | β |
|---|---|---|---|
| N (literal) | (N³ − 1)/3 | (1 − 3U)^{1/3} | **−1/2** |
| ln N | (N² − 1)/2 | (1 − 2U)^{1/2} | 0 |
| −1/N | N − 1 | 1 − U | +1/2 |
| −1/(2N²) (coordinate slowness n = 1/N²) | ln N | e^{−U} | **+1** |

For every f, the linearised speed of δN on a uniform background N₀ is N₀², the coordinate light speed, so all four share the light cone at linear order.

## Parts

**49A. Static mass on the lattice (3D).**
- **Equation and source.** Solve Σ_y N_xN_y [f(N_y) − f(N_x)] = h² s_x on a cube with L = 128 and h = 1. The source is a Gaussian core of radius 1.5 with total S.
- **Boundary.** Dirichlet on the outer faces, F(N) = −(S/4π)/r. This follows from Gauss's law and is not fitted.
- **Solver.** Newton–Krylov (scipy) to a residual below 1e-11.
- **Strengths.** Choose S so that U(r = 6) = 0.02, 0.05, 0.10 and 0.15.
- **Fit.** Over 6 ≤ r ≤ 36, using sites on the three axes and the body diagonal, fit N² = 1 − 2A/r + 2βA²/r² + c₃/r³ with A, β, c₃ free. Use the strengths jointly, with A free per strength and β common.
- **Variables.** All four f.

**49B. One cone (1D chain through a lapse well).**
- **Background.** A chain at impact b = 5 across the static f = N well, N(x) = N₄₉A(√(x² + b²)), at the U(6) = 0.10 strength. Use the 49A radial profile, interpolated.
- **Fields.** Light obeys the 37B rule ∂_τ² u = Σ N_xN_y (u_y − u_x). The firing-rate perturbation δN obeys the linearisation of the f = N rule about the static well. That linearisation includes the gradient terms Σ N_x δN_y (N_y − N_x) + δN_x Σ N_y (N_y − N_x).
- **Packets.** Gaussian packets with k₀ = 0.25, 0.5, 1.0 and width 8/k₀, launched at x = −120. They are evolved in coordinate time with a small fixed dt (RK4), and both are read at x = +120.
- **Readouts.** The arrival-time difference between the light and δN centroids, and each packet's Shapiro delay relative to a flat chain.

**Separate path.** For 49A, a continuum radial solution. Integrate (1/r²)(r² n² f′(n) N′)′ = s(r) by solve_bvp with the same Gaussian source and the far-field condition, apply the same fit window and fit to get β, and also check against the analytic exterior form above.

## Predictions (scored)

- **P1 (lattice β).** The fitted lattice β equals the table value within ±0.05 for each f: −0.5, 0, +0.5, +1.0.
- **P2 (Gauss).** The fitted A equals S/4π within 2% at every strength, for every f.
- **P3 (continuum recompute).** The continuum solution gives the table β within ±0.02 with the same fit, and agrees with the lattice N(r) to within 1e-3 over the fit window.
- **P4 (one cone).** The light and δN Shapiro delays agree within 1% at k₀ = 1. Their arrival-time difference falls with k₀ with a log-slope between −1.5 and −2.5, since the gradient terms act as a weak potential of order ∇²N/N and matter only at long wavelength.

## What the outcomes mean

- **If f = N gives β = −1/2 (P1):** the literal idea 1 is ruled out by perihelion advance. It would give 1.5 times the GR value (Mercury about 64″ per century instead of 43″).
  - The rule still shares the light cone, and it derives a wave equation for N rather than prescribing one.
  - But the quantity a coincidence passes on would have to be the coordinate slowness 1/N² (equivalently, ln N harmonic, the exponential metric) to match GR at this order.
  - The next question would be whether anything in the event picture singles out 1/N² as the transported quantity, for example the refractive index that light itself sees.
- **If the lattice β departs from the continuum value:** the discreteness itself enters at second order. That would be new, and I would report its size against the lattice spacing.

## Outputs

`src/exgc49.py`, `src/exgc49_recompute.py`, `results/EX-GC49/` (A.json, B.json, recompute.json, note). LOG row and HANDOFF section.
