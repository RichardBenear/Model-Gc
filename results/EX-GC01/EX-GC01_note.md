# EX-GC01 (EXPLORATORY): action and reaction for unequal Model Gc particles

Run by Claude, 2026-11-18, under the PI's go-ahead ("yes, run that as an exploration"). About 0.2 core-hours of the
2 core-hour cap. Linux cloud container. Exploratory: this generates hypotheses and supports no claim on its own.

## Result
- **The linear far-field law holds particle by particle.** Each particle approaches at v_i = μ_i M_j/(4πd²) within
  1–2% (0.977–0.991 of prediction), where μ_i is its own mobility and M_j is its partner's consumption. The small
  shortfall is the same one seen in MG195.
- **Weighted action and reaction holds when mobilities match.** We built a partner with 13% more consumption but the
  same mobility within 0.7% (core W diffusion D_W = 0.2 instead of 0.1). The ratio (M_a v_a)/(M_b v_b) is
  1.008–1.009, against 1.007 predicted from the mobility ratio. So the M-weighted centre of the pair stays put to about
  0.1–0.2%: M behaves like a Newtonian mass.
- **It fails, exactly as predicted, when mobilities differ.** We also built a partner with a doubled consumption rate
  k (0.4). Its M is 46% higher and |μ| is 81% higher. The ratio is 0.554–0.556, against 0.552 = μ_a/μ_b predicted.
  The heavier particle moves *more* than the lighter one, and a high-k particle falls 1.8 times faster than a normal
  one in the same field. That is a large violation of the universality of free fall.

## What it means
The Newtonian picture, with M as both the gravitational charge and the inertia, works in Gc only for particles that
share one mobility μ. μ is set by how strongly the core couples to G. Changing the core's size or spread (D_W) leaves
μ almost unchanged. Changing the consumption rate (k) or the W autocatalysis (c) changes μ along with M, and by more
(radial scan, `radial_scan.json`):

| variant | M | μ | μ/M |
|---|---|---|---|
| reference (k 0.2, c 1, D_W 0.1) | 5.318 | −0.02554 | −0.00480 |
| D_W 0.2 | 5.997 | −0.02536 | −0.00423 |
| k 0.4 | 7.752 | −0.04627 | −0.00597 |
| c 1.1 (scan) | 6.723 | −0.03427 | −0.00510 |

Real matter obeys universal free fall to about 1e-15 (MICROSCOPE). For Gc to map onto gravity, all particle types
must share one μ. That is a structural requirement on how particles couple to G, and a new design target: ideally a
symmetry that fixes μ instead of tuning. EX010's 1–2% agreement across +, − and neutral types is the encouraging sign.
The particle types there differ in their charge chemistry, not in their G coupling.

## Method
- **Single particles:** radial 3D solves (the MG195 P1 solver, with per-particle c and D_W) continued from the
  reference profile. All are stable (leading eigenvalues negative apart from the translation mode).
- **Pairs:** a new full-domain axisymmetric solver (`src/exgc01_unequal.py`, `AxiFull`), with no mirror plane.
  - Particle a sits at z = −d/2 and b at +d/2. W-equation parameters are piecewise (a for z < 0, b for z > 0); W is
    effectively zero at the midpoint.
  - Two co-moving speeds come from a bordered Newton solve with two phase conditions. A smooth switch
    w_a = (1 − tanh(z/2))/2 assigns each region's advection to its own particle.
  - Grid as MG195 (h = 0.1 near and between the particles, stretched to 150, G fixed at the outer boundary).
  - The bordered system is solved by block elimination on one sparse LU.
- **Validation.**
  - The equal pair reproduces the MG195 mirror solver: v = 2.18905e-5 against 2.18883e-5 at d = 22 (1e-4), with
    v_a = −v_b to 1e-9.
  - Swapping a and b reproduces the mirrored speeds to 1e-9.
- **Records:** `pairs.json` (all pairs), `radial_*.npz`, `radial_scan.json`. The pair states are kept out of git;
  their checksums are in `pair_states_sha256.txt`. Figure: `paper/figures/EX-GC01_action_reaction.png`.

## Limits
- Speeds come from quasi-static co-moving solves (separation fixed). That is valid because the speeds are about 1e-5.
- One reference particle, three separations, two partners. The switch width is fixed at 2; it acts only where W ≈ 0.
- The motion is still overdamped (velocity, not acceleration). See docs/PHYSICS_MAPPING.md §2.
