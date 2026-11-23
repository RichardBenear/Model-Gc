# EX-GC34E — with a shift link variable, the shared structure reproduces linear GR (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "you can close 34D and go on to 34E". Plan
`spec/EX-GC34E_plan.md` (ab387ca, pushed before the analysis). Code `src/exgc34e.py` → `analysis.json`;
`src/exgc34e_check.py` → `check.json`. Post hoc: `src/exgc34e_posthoc_modes.py` → `posthoc_modes.json`,
`src/exgc34e_posthoc_cutoff.py` → `posthoc_cutoff.json`. Our construction, not LaViolette's Model G.
**PI CLOSED 2026-11-23** (Richard: "go ahead and close 34A to 34E"). **Corrected 2026-11-23** after closure (see
Correction): the vacuum spectrum has a growing mode.

## Correction (2026-11-23, after ChatGPT's review)

**The conclusion below first said the lattice has exactly two wave modes and no extra mode. That was wrong.** The
vacuum equations also have a pair of modes with imaginary frequency, one of which grows exponentially. They were in
the committed `analysis.json` (E3 `all_finite`, directions 110 and 111), but my code kept only real roots, and my post
hoc scan searched real frequencies only. ChatGPT found them with a separate reduction; I confirmed them with
`src/exgc34e_review_unstable.py` → `review_unstable.json` (post hoc, not a re-scoring; E3 P1 stays FAIL as
registered):

| kh | 110 | 111 | generic (0.8, 0.5, 0.33) | near axis (1, 0.05, 0) |
|---:|---:|---:|---:|---:|
| 0.025 | ±2.122i | ±1.415i | ±1.421i | ±1.573i |
| 0.05 | ±2.124i | ±1.416i | ±1.422i | ±1.574i |
| 0.1 | ±2.132i | ±1.420i | ±1.426i | ±1.577i |
| 0.2 | ±2.163i | ±1.436i | ±1.441i | ±1.590i |

Entries are ω/k (convention e^{−iωt}, so +i grows). Each root satisfies all 11 equations (relative residual ≤ 1e-15
for 111 and generic, ≤ 5e-7 where the pencil also has an exact null direction), and the gauge-fixed kernel's smallest
singular value falls to ~1e-18 there from ~1e-11 nearby on the imaginary axis.

- **The growth rate does not shrink as the lattice is refined** (Im ω ≈ 1.4–2.1 k for every kh tested). It is not a
  small lattice correction.
- **It is the re-slicing direction.** The mode has lapse and lengthwise shift only, with no edge change in this gauge
  (|n|/|β∥| = Im ω/k). About 99.98% of its metric perturbation at kh = 0.05 is a continuum coordinate change; the
  curvature-carrying part is ≈ 0.06 (kh)² of it. Where the re-slicing symmetry is exact this direction is pure gauge;
  where it is broken (off the lattice axes) the lattice turns it into a growing solution. Because it grows
  exponentially, its small curvature part would still matter in a time-stepped run.
- **Along the 100 axis this method cannot decide.** The gauge-fixed pencil has three exact null directions there.
  Wavevectors in a lattice plane (110, and (1, 0.05, 0)) have one exact null direction; generic directions have none.
- The two oscillating modes stand: their speeds differ from the 34A light speed (0.999896 at kh = 0.05, the same in
  all three directions) by at most 1.4e-4. ChatGPT quoted 2.5e-4 for the faster 110 branch; I could not reproduce
  that number and have not resolved the difference.

**Revised verdict.** The forced moving-source result (E4: gravitomagnetism and the transverse momentum constraint)
stands. The vacuum spectrum does not match linear GR: besides the two light-speed modes there is an unstable branch
in the re-slicing direction. So 34E does **not** establish GR's constraint structure at finite spacing, and "linear GR
on the lattice" is too strong. A stability study (34F, not yet approved) must come before any wave-propagation or
non-linear work.

The original conclusion is kept below, unchanged, for the record. Its wave and constraint-structure statements are
superseded by this correction.

## Original conclusion (superseded in part, see above)

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

- **Established (linear, flat background):** GR's static and stationary (moving-source, O(v)) linear solutions, with
  lattice corrections of O((kh)²), and two light-speed radiative modes. *(Corrected: the original line also claimed
  GR's constraint structure and only GR's radiative solutions. The growing re-slicing mode in the Correction section
  contradicts both.)*
- **Not established:** non-linear gravity (for example perihelion precession or strong fields), a lattice time step
  (here time is continuous), and whether the weakly broken re-slicing symmetry is harmless once the equations are
  stepped in time and the gauge must be fixed explicitly. Regge-type lattice gravity is known to have this issue at
  the non-linear level.
- **Model link:** none of 34A–34E uses the Gc reaction–diffusion particles. The knots (GC26) are the matter here. How
  Model Gc's conserved-G particle would sit on this structure is open.

## Compute

About 0.2 core-hour of the 2 cap: registered analysis 5 CPU-s, check seconds, post hoc mode scans about 10 minutes.
Linux cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1, mpmath 1.3.0.

## Suggested next (not approved; superseded by the Correction: 34F should be a stability study)

One option is to step back and consolidate. Write 34A–34E into a white-paper v0.5 section ("one shared structure:
universal coupling, light bending, gravitomagnetism, waves"), with ChatGPT doing the recomputation and review the GC34
outline foresaw. Or, scientifically, **34F**: time-step the linear lattice equations (leapfrog), check that a
gravitational-wave packet propagates stably at light speed, and see whether the weakly fixed re-slicing direction
grows. That is the first place the non-exact symmetry could bite.
