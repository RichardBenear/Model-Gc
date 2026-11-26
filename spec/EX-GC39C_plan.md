# EX-GC39C — does a free lapse field radiate? (scalar radiation from an orbiting pair)

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread (2026-11-26 04:53 UTC). Written before any run. Cap 4 core-h. Design: [EX-GC39_design.md](EX-GC39_design.md). This is our construction, not LaViolette's Model G.

## Question
If N is a free relativistic field (37C), with n = N − 1 obeying n_tt − ∇²n = −4πGρ, so that the static solution is N = 1 − GM/r as in 37B, how much power does a binary radiate into n? How does that compare with GR's quadrupole formula, which binary pulsars confirm to about 0.2%?

## Prediction (analytic, recorded before the run)
- Lagrangian (1/4πG)·½(ṅ² − |∇n|²) − ρn. This gives attraction and positive wave energy, with flux S = −(1/4πG)ṅ∇n.
- Far field: n = −(G/r)[M + ½n̂_in̂_jQ̈_ij(t − r) + …]. No dipole term, because the COM is fixed and the coupling is universal.
- **Power:** P_s = (G/16π)∮(n̂n̂Q⃛)² dΩ = (G/30)Q⃛_ijQ⃛_ij, since the trace of Q is constant on a circular orbit.
- **Equal-mass circular binary** (masses M, radius a, frequency Ω): P_s = (64/15)G M² a⁴ Ω⁶. GR: P_GR = (128/5)G M² a⁴ Ω⁶. So **P_s/P_GR = 1/6**.
- **Exact (non-multipole) value for a prescribed source:** P = (Gω²/8π)∮|ρ̃_ω(ω n̂)|²dΩ, summed over harmonics. This is used as the comparator for finite source size.

## Lattice run
- 3D cubic lattice, h = 1, 160³, leapfrog dt = 0.25, absorbing sponge in the outer 20 cells.
- Prescribed sources: Gaussians of width σ = 1.0, total mass M each, on a circular orbit of radius a = 2.5 about the box centre, Ω = 0.08 (radiated wavelength 2π/(2Ω) ≈ 39).
- Unequal-mass case: M1 = 2M2 with radii in the inverse ratio (COM fixed), to check that there is no dipole.
- Readout: the time-averaged flux through spheres of radius 30, 40 and 50, over the last two periods after transients (t from 150 to 150 + 4π/Ω).

## Predictions (scored)
- **P1:** the lattice power matches the exact comparator within 5% (lattice dispersion at k = 0.16 is below 0.3%). This checks the solver.
- **P2:** the unequal-mass power matches the quadrupole formula with Q from both masses within 5%: no dipole radiation.
- **P3 (the physics):** P_s/P_GR = 1/6 within 5% (with the finite-size factor). Consequences, stated before the run:
  - **If gravity were only N (no tensor modes),** orbits would decay at 1/6 of the observed rate, which is ruled out.
  - **If N coexists with GR's tensor waves,** there is an extra 17% decay, ruled out by about 80× the binary-pulsar precision.
  - **Either way,** N cannot be a free radiating field with full gravitational coupling. It must be constrained, as in H-A, or its waves must not couple at this strength.

## Checks
- **G1:** static source: the lattice solution matches −GM/r for r = 10–40 within 1%.
- **G2:** flux independent of sphere radius within 3% (energy conservation outside the source).

## Outputs
`src/exgc39c.py`, `results/EX-GC39C/`.
