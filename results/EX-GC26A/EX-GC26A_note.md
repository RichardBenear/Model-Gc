# EX-GC26A — a conservative 3D topological soliton (EXPLORATORY)

Claude, 2026-11-20. PI go-ahead: “Go ahead with your proposed 26A test.” Plan and pass criteria logged in
`explore/LOG.md` before the production runs (commit cbaa488). Code `src/exgc26a_skyrmion.py`, figure
`src/exgc26a_figure.py` → `paper/figures/exgc26a_skyrmion.png`, records `results/EX-GC26A/*.json`. The model is the
plain Skyrme model (Skyrme 1961), not Model Gc: no chemistry and no G. The discretisation, tests and criteria are our
construction.

## Conclusion

**The physics works; three of the registered numerical criteria did not.** A degree-1 Skyrmion on the lattice is a
stable, lossless particle. Its rest energy converges toward the independent radial value, its total energy is
conserved to 3 parts in 10¹⁰ while it is shaken, its knot number holds, and after a kick it settles back to the same
rest energy (to 2 parts in 10⁷). A no-knot start of similar size simply dissolves. But the registered criteria were
met for 5 of 8: the lattice degree at h = 0.2, the Derrick energy balance, and the lopsided-start test fell short, and
total momentum drifted slightly more than allowed. All three numerical shortfalls shrink with finer grids or more
iterations as far as checked; none is a sign that the particle is unstable. The lattice itself breaks the frame
independence a little (momentum drift), which matters for 26B.

## Results against the registered criteria

| Criterion | Registered requirement | Result | Verdict |
|---|---|---|---|
| C1 degree | \|B\| within 0.05 of 1 at h = 0.2 and 0.15 | 0.929 (h = 0.2), 0.962 (h = 0.15) | **Not met** at h = 0.2; met at 0.15 |
| C1 energy | finite, positive; h² extrapolation within 1.5% of radial; box ±7 vs ±5 within 1% | 1.2061 (0.2), 1.2212 (0.15) → 1.2406, +0.75% from radial 1.2314; box ±7 1.2024 (−0.31%) | Met |
| C2 Derrick balance | E₂/E₄ within 2% of 1 at h = 0.15 | 0.911 | **Not met** |
| C3 lopsided start | same E within 0.5%, same degree | E 1.3105 (+8.7%) after the 20 000-iteration limit; degree 0.931 | **Not met** (not converged) |
| C4 no-knot control | E < 1e-3 | 6e-34 | Met |
| C5 energy ledger | conserved within 1e-3 to T = 6 | 2.6e-10 | Met |
| C6 momentum ledger | \|ΔP\| < 1e-3 E₀ for t ≤ 3 | max 1.31e-3 E₀ (t = 2); 3.8e-3 by t = 6 | **Not met** |
| C7 degree in the kick | \|ΔB\| < 0.02 | 0.0053 | Met |
| C8 recovery | back to E₀ within 0.5%, same degree | +2.3e-7, degree unchanged | Met |

Energies are E/12π² in Skyrme units, where the topological bound is E ≥ \|B\|. The kick added kinetic energy equal to
5% of E₀ with a random smooth velocity field (net momentum 0.126 E₀).

## Reading the shortfalls

- **The lattice degree is not an integer.** On a grid the computed degree is a sum that approaches 1 as the grid is
  refined: 0.875, 0.929, 0.962 at h = 0.25, 0.2, 0.15, shrinking about as h². The 0.05 tolerance was too tight for
  h = 0.2. The degree did not drift during the dynamics (C7).
- **The Derrick balance is off by lattice and box effects, not by poor convergence.** A post-hoc check continued the
  h = 0.2 relaxation for 6000 more iterations: the energy changed by 1.4e-9 and E₂/E₄ stayed at 0.861. The ratio rises
  with finer grids (0.776, 0.861, 0.911) and with a bigger box (0.879 at ±7), so it is heading toward 1 but slowly.
- **The lopsided start had not finished relaxing.** It began with a lot of noise (E = 9.3) and ended at 1.3105 with the
  right degree when it hit the iteration limit. Whether it would reach 1.2061 is not shown here. Not rerun, because the
  compute cap was reached.
- **Momentum drifts by about 1% of the kick's momentum.** A grid is not translation invariant, so the field exchanges a
  little momentum with it. Energy is still conserved to 1e-10. This is the lattice's own rest frame (as in EX-GC25) and
  needs to be measured as a function of h in 26B before inertia is measured.

## Other findings

- **Plain central differences fail.** In solver tests the knot escaped through the grid's checkerboard mode, which costs
  no energy with central differences. Averaging all-forward and all-backward differences removes that mode.
- **The knot needs about five grid points across.** In solver tests the knot collapsed through the grid at h = 0.3 and
  0.4: there the discrete energy keeps falling as it shrinks. At h = 0.25 it survives (degree 0.875, E 1.184).
- **Separate path.** The radial shooting calculation gives E/12π² = 1.23136 (literature 1.2315) and an exact Derrick
  balance (E₂ fraction 0.49997). It shares no code with the 3D solver.

## Compute

Production runs 2.02 core-hours (`cpu_s` in the records), against the 2 core-hour cap: a small overrun, mainly
because the control and lopsided runs ran to the iteration limit. Solver tests before the log entry and the post-hoc
convergence check added about 0.4 core-hours. Total about 2.4 core-hours in the cloud, over the cap; recorded here as
a deviation. Platform: Linux cloud container, Python 3, numpy 2.x.

## For 26B (PI decision)
The particle is ready for boosts if the grid is h ≤ 0.15. Before measuring inertia, 26B should measure the lattice's
momentum leakage as a function of h, since it acts like a small drag. 26B will need more than the 2 core-hour cloud
cap; the PI's PC is the better place.

## Correction (2026-11-20, after ChatGPT's review)
The scoring "5 of 8" above counted C1 as met on its energy part. As registered, C1 also required the degree condition
at h = 0.2, which failed, so C1 is a **partial pass** and **4 of 8 criteria were fully met** (C4, C5, C7, C8). The
results are unchanged. EX-GC26A2 later showed with a geometric degree that the particle's degree is exactly −1 in every
resolved state; the registered C1 verdict stays as scored.
