# EX-GC26B — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-21. Reviewed commits d4c98a1 (plan) and 0c0004c (results), by ChatGPT. Plain massless
Skyrme model; the numerical construction is ours.

## Verdict

Sound and honestly reported. I agree that the full 1% precision gate is **not met** and that precision 26E should
stay blocked. But I read the cause more narrowly than the note does. The boosted knot carries the right energy and the
right momentum. What goes wrong in the axial runs is that its energy does not flow at the speed its momentum implies.
That is a known property of the two-orientation stencil, not a failure of relativistic inertia. A targeted 26B2 with a
mirror-symmetric stencil and a longer coast should be able to fix it.

## What I checked

1. **Recomputation of all 15 runs from the committed time series**, with my own fitting code (not `exgc26b_figure.py`):
   measured speed, M/E₀, H/(γE₀) and transverse momentum match the summaries to every printed digit. The dt-halved run
   matches the dt = 0.05h run to 10 digits, as stated.
2. **Independent rerun with my original EX-GC26A NumPy solver** (no C++ kernel, no A2 code), from my own 26A state
   (h = 0.2, box ±5, not ±7). I built φ_t by differencing the Lorentz-transformed field in time, not by differentiating
   the spline in space. u = 0.3, T = 0.8. Script and outputs: `review/EX-GC26B_check/`.

   | direction | quantity | mine (L = 5) | ChatGPT (L = 7) |
   |---|---|---:|---:|
   | axis | M/E₀ | 1.0222 | 1.0215 |
   | axis | centroid speed ÷ (P/H) | 0.9817 | 0.9822 |
   | axis | transverse P/E₀ | 3.8e-3 | 3.7e-3 |
   | diagonal | M/E₀ | 1.0044 | 1.0043 |
   | diagonal | centroid speed ÷ (P/H) | 0.9988 | 0.9988 |

   The axial excess and the axial/diagonal split reproduce in a separate code and a smaller box. Same solver author
   for the 26A code, so this confirms the numbers, not the stencil's physics.

## The key diagnostic: does energy move at P/H?

In the continuum, the energy centroid moves at exactly P/H (the rate of change of ∫x T⁰⁰ is the momentum). So the
ratio of the measured centroid speed to P/H tests the lattice directly, independent of E₀ and of the launch:

| h | u | axis: speed ÷ (P/H) | diagonal: speed ÷ (P/H) |
|---:|---:|---:|---:|
| 0.2 | 0.1 / 0.3 / 0.5 | 0.986 / 0.982 / 0.971 | 1.000 / 0.999 / 0.995 |
| 0.15 | 0.1 / 0.3 / 0.5 | 0.993 / 0.991 / 0.985 | 1.000 / 0.999 / 0.998 |

Meanwhile the momentum and energy themselves are close to the Lorentz values in every fine-grid case:
P/(γuE₀) = 0.9995 to 1.0085 and H/(γE₀) = 1.0001 to 1.0042. So:

- The mass "error" is almost entirely a low measured speed, and the low speed is present from the first frames
  (first-half and second-half fits agree), not a slow build-up.
- It is large only for axial boosts, which are also the only runs with transverse momentum. The averaged
  all-forward/all-backward stencil is symmetric under swapping axes and under full inversion, but not under a single
  mirror (y → −y). That allows sideways momentum in an axial boost and lets the lattice energy current differ from the
  canonical momentum along the axes. Diagonal boosts are protected by the swap symmetry, and there the identity holds
  to 0.01 to 0.5%.
- The axial deviation shrinks by about 2× from h = 0.2 to 0.15, roughly as h².

**Two-point h² extrapolation, indication only.** ChatGPT rightly declined to extrapolate from two resolutions. As a
rough indication: M/E₀ extrapolates to 0.993–0.998 (axis) and 0.998–1.000 (diagonal) across the three speeds.

## Other comments

- **Momentum loss here is the lattice, not the box.** The flux correction barely changes the regional residual. The
  h = 0.2 axial 0.5c regional residual equals the global ΔP (both 4.83e-3 E₀), and L = 9 changes almost nothing. This
  revises my EX-GC26A2 comment for fast knots: the box dominated the slow, kicked case, but at 0.3 to 0.5c the
  lattice's lack of exact translation symmetry dominates. The centroid also decelerates slightly (h = 0.2 axial 0.5c:
  about −0.0065; diagonal −0.0013). This is the lattice acting as a rest frame, as anticipated in EX-GC25. At T = 0.8
  (the knot moves under half its own width) the size of that drag is not yet measured.
- **Static tensor caveat is well taken.** The centred and matched estimators disagree by up to 4% at h = 0.2. That
  shows the "mass" depends on how translation is represented on the lattice. The dynamic P/H test above avoids that
  choice.
- **Deferring the pion mass is reasonable** for comparability with 26A/A2.
- **Protocol:** clean. The plan reached main 16 minutes before the results. Production used 0.56 core-hours of the 2
  cap. My check used about 0.1 core-hours on the cloud container.

## Recommendation for 26B2

Agree with ChatGPT's targeted 26B2, with these specifics:

1. A reflection-symmetric stencil, averaging all eight one-sided orientations, with its own consistent kinetic term.
   Re-validate it: energy conservation, geometric degree, and a checkerboard check.
2. Use **centroid speed ÷ (P/H) within 0.5%, and transverse P/E₀ < 1e-3, along the axis** as the primary
   stencil-quality gate. It is the cleanest discriminator found so far.
3. A longer coast (T ≥ 5 to 10) at u ≤ 0.3, to measure the lattice drag rate directly, since 26E needs the knot to
   move freely over many widths.
4. A force-response test (constant weak applied force, measure acceleration), as ChatGPT suggests. For 26E this is
   the quantity that matters, and gravity tests are low speed.
5. Optionally, a third resolution (h = 0.125) for one axial case to confirm the h² trend.

For now, the fine-grid low-speed operating point (h = 0.15, u ≤ 0.1: axial +0.73%, diagonal +0.02%) is good enough
for qualitative exploratory work, as the note says. ChatGPT should own 26B2 (it has the kernel), with me reviewing,
unless the PI prefers otherwise.
