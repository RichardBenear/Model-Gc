# EX-GC34E — with a shift link variable, the shared structure reproduces linear GR (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "you can close 34D and go on to 34E". Plan
`spec/EX-GC34E_plan.md` (ab387ca, pushed before the analysis). Code `src/exgc34e.py` → `analysis.json`;
`src/exgc34e_check.py` → `check.json`. Post hoc: `src/exgc34e_posthoc_modes.py` → `posthoc_modes.json`,
`src/exgc34e_posthoc_cutoff.py` → `posthoc_cutoff.json`. Our construction, not LaViolette's Model G.

## Conclusion

**Adding one more link variable, a shift (how fast the lattice vertices move relative to the clock-normal
observers), turns the 34D structure into linear general relativity on the lattice.** The kinetic energy is the
lattice version of GR's, built from each tetrahedron's strain rate and divided by the lapse as for matter.

- **Moving sources:** the field of a moving mass, including its gravitomagnetic part, matches linearised GR to
  4e-4 at kh = 0.05. The error shrinks as k², so it is a lattice-spacing effect. The case 34D's no-shift structure
  got wrong (motion perpendicular to the wavevector) now satisfies the momentum constraint to 1e-15.
- **Gravitational waves:** there are exactly two wave modes, and they travel at the same speed as light on the same
  lattice, within 1.4e-4 at kh = 0.05. No extra scalar mode was found anywhere from very low frequency up to the
  lattice cutoff. The two polarisations split slightly off the lattice axes (2e-4 at kh = 0.05, growing as k²).
  That is lattice anisotropy.
- **Re-slicing time is an exact lattice symmetry only for waves along a lattice axis.** In other directions it is
  broken at O((kh)²) relative to the physical modes. In this linear analysis the breaking produces no extra light
  mode. It does leave one field direction that is only weakly fixed, which could matter in a time-stepped or
  non-linear version.

So the chain 34A → 34E gives this result for our construction: one shared structure (lapse, shift, Regge edge
lengths) that every field moves through, whose own energy follows the same clock rule as matter's. That reproduces
universal free fall, the Newtonian limit, GR light bending, gravitomagnetism, and light-speed gravitational waves with
two polarisations. All of this is at linear order on a flat background.

## Registered results

| part | criterion | result | verdict |
|---|---|---|---|
| E1 | displacement modes null: ‖MG‖/‖M‖ ≤ 1e-10 | 3.7e-16 | PASS |
| E1 | static limit = 34D kernel ≤ 1e-12 | 0 | PASS |
| E1 | real-space kinetic sum vs Fourier form ≤ 1e-10 | 1.4e-15 (8³ lattice) | PASS |
| E2 | re-slicing symmetry (reported) | exact along axes (1e-17); off-axis, smallest singular value ∝ (kh)⁴ of the largest, i.e. ∝ (kh)² of the physical scale | as expected (approximate) |
| E3 P1 | exactly two modes with ωh < 0.5 at kh = 0.05, 0.1 along 100, 110, 111 | found 0, 1, 2 | **FAIL as registered** (method, see below) |
| E3 P2 | \|ω/k − 1\| ≤ 1e-3 and \|c_gw − c_light\| ≤ 1e-3 at kh = 0.05 | 110: 0.99984 (c_light 0.99990); 111: 1.00002; 100: none found | PASS where measured |
| E4 P4 | O(v) Riemann vs GR ≤ 1e-2 at kh = 0.05, shrinking with k (9 cases) | max 4.2e-4; 1.0e-4 → 4.2e-4 → 1.7e-3 at kh = 0.025, 0.05, 0.1 | PASS |
| E4 P5 | v ⊥ k: finite, momentum constraint ≤ 1e-10 | residual ≤ 1.2e-15; Riemann as P4 | PASS |

**The registered wave-mode method failed, not the physics.** I found the modes by turning the quadratic eigenproblem
into a matrix pencil after fixing the displacement gauge. That only works if the pencil is non-singular. Where the
lattice keeps the re-slicing symmetry exactly (waves along an axis), the pencil is singular, and the eigensolver
returns no reliable finite modes. Off-axis the symmetry is nearly exact, which loses modes partially. The plan should
have removed that direction too. P1 stays FAIL as registered.

**Post hoc mode search (not a re-scoring).** I projected out the displacement gauge at each frequency and scanned for
the frequencies where a further singular value of the kernel drops to zero (< 1e-11 relative, or to the re-slicing
floor off-axis):

| kh | 100 | 110 | 111 | light (34A field) |
|---:|---|---|---|---:|
| 0.05 | 0.999757, 0.999948 | 0.999839, 1.000029 | 1.000017 (pair) | 0.999896 |
| 0.1 | 0.999029, 0.999792 | 0.999355, 1.000114 | 1.000070 (pair) | 0.999583 |
| 0.2 | 0.996134, 0.999164 | 0.997428, 1.000453 | 1.000284 (pair) | 0.998334 |

Speeds are ω/k. Two modes in every direction, and nothing else. A log scan of the weakly fixed (re-slicing)
direction found no zero crossing from ω = 1e-4 k up to ωh = 2.5. Its singular value stays smooth, within a factor of
about 1 of its static value up to ω = 3k and rising above that.

## What's established, and what isn't

- **Established (linear, flat background):** the shared structure with lapse, shift and Regge edges has GR's
  constraint structure. Lapse and shift are multipliers (Hamiltonian and momentum constraints). It has GR's static,
  stationary and radiative linear solutions, with lattice corrections of O((kh)²).
- **Not established:** non-linear gravity (for example perihelion precession or strong fields), a lattice time step
  (here time is continuous), and whether the weakly broken re-slicing symmetry is harmless once the equations are
  stepped in time and the gauge must be fixed explicitly. Regge-type lattice gravity is known to have this issue at
  the non-linear level.
- **Model link:** none of 34A–34E uses the Gc reaction–diffusion particles. The knots (GC26) are the matter here. How
  Model Gc's conserved-G particle would sit on this structure is open.

## Compute

About 0.2 core-hour of the 2 cap: registered analysis 5 CPU-s, check seconds, post hoc mode scans about 10 minutes.
Linux cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1, mpmath 1.3.0.

## Suggested next (not approved)

One option is to step back and consolidate. Write 34A–34E into a white-paper v0.5 section ("one shared structure:
universal coupling, light bending, gravitomagnetism, waves"), with ChatGPT doing the recomputation and review the GC34
outline foresaw. Or, scientifically, **34F**: time-step the linear lattice equations (leapfrog), check that a
gravitational-wave packet propagates stably at light speed, and see whether the weakly fixed re-slicing direction
grows. That is the first place the non-exact symmetry could bite.
