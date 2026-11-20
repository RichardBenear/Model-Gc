# EX-GC25 — Short look at a medium with no rest frame (EXPLORATORY)

Claude, 2026-11-20. PI go-ahead: “The short look at b) is my preference. Go ahead and run this.” Logged in
`explore/LOG.md` before the run (34ae9b8). Code `src/exgc25_framefree.py`, record `results/EX-GC25/EX-GC25_runs.json`,
figure `paper/figures/exgc25_framefree.png`. 13 CPU-seconds. Model Gc and everything here is our construction.

## Conclusion

Gravity can survive without a rest frame, and it comes out better than in Model Gc: inverse-square, attractive,
and with universal free fall built in. But it is no longer gravity from consumption. In a frame-free medium, G
cannot be a conserved carrier that particles absorb. It has to be a field that particles disturb in proportion to
their energy. Inertia then comes for free, and so does the equivalence principle. This is a different model, not a
repair of Model Gc.

## The argument

1. **Why consumption forces a rest frame.** A conserved G has a density and a flux that together form one
   four-vector. A uniform background density G∞ is that vector pointing purely in time, which is true in one frame
   only. A particle moving through it at speed u meets a flux G∞u, and absorbing that flux costs it momentum. That is
   exactly the head-on drag b = S found in EX-GC17 and EX-GC24A. So with a conserved, absorbed carrier, the drag can't
   be removed by any choice of medium dynamics. This is general, not specific to our closure.
2. **What can have no rest frame.** A constant scalar field value looks the same in every frame. So let G be a scalar
   field with vacuum value G∞ and a wave equation (speed c, no friction). A particle at rest disturbs it, and the
   static disturbance obeys Laplace's equation, so the 1/r deficit and the 1/d² force survive. For a scalar field,
   like sources attract. A particle moving at constant speed carries a Lorentz-contracted copy of its field, so
   there is nothing to drag on.
3. **What the field must couple to.** If G couples to the particle's energy (technically to the trace of its
   energy–momentum tensor, Nordström's 1913 scalar theory), then a theorem by von Laue says that for any stable,
   localised, static object the pressure terms integrate to zero. The gravitational charge is then exactly the rest
   energy, which is also the inertial mass. This holds for every particle type, and for bound composites including
   their binding energy. So the failures of EX-GC01 to EX-GC04 (non-universal fall across types and composites) are
   fixed automatically. If instead G couples to anything other than energy (like the absorption rate S in Model Gc),
   the fall rate is Q/M and varies by type.

## The 1+1 D check

Relativistic φ⁴ kinks stand in for particles (they are not Gc particles). Four kinks, A to D, with rest energies
spread 8× and widths spread 2.7×, start at rest in a weak uniform potential Φ = g x.

| Coupling | A | B | C | D | Prediction |
|---|---|---|---|---|---|
| to energy (Nordström) | 0.9995 | 0.9987 | 0.9987 | 0.9995 | 1 for all |
| to a non-energy charge | 1.0019 | 0.2499 | 0.0625 | 0.7841 | 1, 0.25, 0.0625, 0.790 |

The values are fall acceleration divided by g. Coupled to energy, all four fall together to 0.13%. Coupled to
anything else, they spread by a factor 16, exactly as Q/M predicts (within 0.8%).

**Coasting.** A kink launched at half the wave speed, with no potential, keeps its speed. Over 200 time units it
loses 2e-5 to 2e-4 of its speed, shrinking as the lattice gets finer. Compare a coasting time of about 0.04 for the
Model Gc particle in EX-GC24A. The exception is a kink only 1.8 lattice spacings wide, which loses 14%: the lattice
itself is a rest frame, and a particle about as small as the lattice spacing feels it.

## Costs and open points

- **G is no longer conserved or consumed.** Model Gc's physical picture (a carrier drawn into the core) goes. The
  1/r field is a static disturbance, not a flow. The name Gc would no longer fit.
- **Particles must be conservative.** Friction picks a rest frame, so the dissipative reaction–diffusion particles
  (Brusselator core, Turing-type) can't live in this medium. The particle has to be a conservative soliton. The
  lattice knot (EX-GC12 to 15) is the natural candidate; in 3D it would need a Skyrme-type term to be stable.
- **Scalar gravity is only part of gravity.** Nordström's theory gives Newton's law and the equivalence principle,
  but it is known to fail later tests: it predicts no bending of light and the wrong perihelion shift for Mercury
  (a backward shift, −1/6 of the observed one). A full theory needs a tensor (spin-2) part.
- **Discreteness is a rest frame.** On a lattice or hypergraph, drag returns for particles near the grid size. The
  1+1 D result suggests it is negligible once the particle spans more than about 4 spacings. Causal invariance is the
  property claimed to restore frame independence at large scales in Wolfram-style models.
- **The 1+1 D check is limited.** It tests a test body in a fixed potential. It does not test the particle sourcing
  the field, 3D, or two particles attracting.

## Possible next step (needs the PI's decision)

EX-GC26: a 3D relativistic lattice knot with a dynamic scalar field G sourced by the knot's energy. Check (i) 1/d²
attraction between two knots, (ii) universal fall for two knot sizes, (iii) coasting with no drag. The alternative
is to keep Model Gc's consumption picture and pursue internal inertia (option a), accepting the drag.
