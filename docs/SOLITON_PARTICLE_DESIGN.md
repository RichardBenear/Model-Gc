# Exchange sign via a soliton particle (design note, 2026-11-19)

> **Later outcome (2026-11-22):** the soliton route was carried out in 3D as EX-GC26 (Skyrme model with Nordström
> gravity); see Section 8 of `paper/Model_Gc_white_paper_v0.4.md`. This note is kept unchanged below.

Author: Claude, at the PI's request (decision card, 2026-11-19: "Particle as a field soliton"). Design only. Points
are marked **established** (standard mathematics or our records), **design** (our choice) or **speculative**. Our
construction, not LaViolette's.

## Why a soliton (established)
- EX-GC11 Step 0: in EX-GC10's frame-field tether, exchanging two particles cannot give a sign. The field carries
  no belt identity, and a twisted core is pinned by the no-cut rule.
- The field version of the spin–statistics link (Finkelstein & Rubinstein, J. Math. Phys. 9, 1762 (1968); Williams
  1970) needs the particle to be a texture of the field itself: a map from space (with the far field uniform) to
  SU(2) of degree B (π₃(S³) = Z). In each degree sector, loops of configurations have a Z₂ (π₄(S³) = Z₂). For odd B:
  - a 2π rotation of the soliton is a non-contractible loop (the 720° property);
  - exchanging two such solitons is homotopic to a 2π rotation of one of them (the fermion sign).

  For even B, both loops are trivial.
- Skyrme's nucleon is this kind of object.

## The discrete version (design)
- Keep EX-GC08/10's lattice: one frame (a cube rotation) per node, with the lift rule (steps 90° and 120°, short
  SU(2) lifts must compose). Because neighbouring lifts are always short, the field has a consistent SU(2) lift q(x)
  in the 48-element binary octahedral group, fixed to +1 at the outer face.
- **Degree.** Split each lattice cube into 6 tetrahedra (Kuhn). Each tetrahedron maps to the spherical tetrahedron
  spanned by its four lifts. Count, with signs, the tetrahedra that cover a generic point of S³.
- **Soliton.** The hedgehog q = cos F(r) + sin F(r) r̂·(i, j, k), with F falling from π at the centre to 0 at
  r = R, snapped to the nearest of the 48 lifts. A prototype check (R = 8, 10, 11 on 25³) gives:
  - neighbouring lifts are always within 120°, so the texture obeys the step rule;
  - the degree comes out as ±1 for random test points.
- **Why it differs from EX-GC10.** EX-GC10's twisted state wraps a circle (degree 0), so it is "object plus belt".
  The hedgehog wraps all of S³, so the particle is the knot itself, with no separate object.

## Stages (proposed)
1. **Stability (cheap, frames only).** Does the lattice dynamics conserve the degree? The test runs Metropolis on
   mismatched edges at fixed T_f with the lift rule.
   - Prediction: a degree-1 texture cannot reach zero energy.
   - Control: the pairwise rule, which should let it unwind.
   - Risk: on a lattice a texture can collapse through one cell, which is the discrete version of Derrick's
     instability. If that happens, the degree is not protected, and a size-setting ingredient would be needed.
2. **The loop test (the actual sign).** Build the two loops as explicit sequences of allowed moves: a 2π rotation of
   one texture, and an exchange of two textures. Compute their Z₂ class. The class is the π₄(S³) invariant of the
   loop, read off from the framed preimage (Pontryagin–Thom) of a regular point in S¹ × S³. This is the hard part.
   It needs a careful discrete construction, so it will get its own entry after Stage 1.
3. **Tie to Gc (speculative).** Let the texture's core be the G consumer, so the soliton is also the gravitating
   particle. One option: G consumption is switched on where the lift is near −1, i.e. at the texture centre. Then
   check that the far field stays 1/r.

## Caveats (speculative)
- Route 2 was set aside because it is unclear what makes a twisted field in nature. Here the field is the graph's
  own edge ordering (the rotation system), as in route 3, so no extra substance is added. The open question moves to
  what creates a knotted rotation system in the first place.
- Even if the loop classes come out right, a classical Z₂ of loops becomes a quantum sign only once the programme
  has a quantum layer (amplitudes over histories). These stages would show the topology is present, not that the
  particles are fermions.
