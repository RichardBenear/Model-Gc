# ChatGPT review of EX-GC34F and EX-GC34G

2026-11-23 | Independent review of the exploratory results

**Verdict.** GC34F confirms the instability of the original GC34E vacuum equations. GC34G found no repair among the seven registered candidates that met all requirements. The static and moving-source benchmarks survive. A stable lattice model with just the two expected long-wavelength gravitational-wave modes remains open.

## GC34F

**Time evolution.** The disturbance grows at the predicted rate in two tested directions, with the lapse and momentum equations satisfied and growth visible in reconstructed tidal curvature. I recalculated the exact implicit-midpoint growth rate, `log((1 + γ dt/2)/(1 − γ dt/2))/dt`, for the four recorded direction/timestep combinations. It matches the reported amplitude growth to within 4 × 10⁻¹⁰ in the run’s time units. This checks the evolution of the predicted mode in the same equations. In the longer post hoc run, the prescribed seed overtook the gravitational wave’s tidal signal after about three wave periods in the generic direction. The registered window ended before crossover.

**Refinement.** At fixed physical wavelength, the sampled growth rates approach a nonzero value as spacing decreases, while the mode becomes closer to a continuum coordinate change. The observable contamination time depends on the physical curvature of the initial seed and the duration of evolution. Restrict the conclusion to the tested discretization and distinguish the deliberately seeded run from claims about excitation by round-off or instability of the continuum Einstein equations.

**Mechanism.** The results implicate the balance between the discrete kinetic and curvature terms in a direction that is a coordinate freedom in continuum GR. The proposed minimal model misses the generic-direction growth rate by about 35%, so the registered mechanism test failed. The kinetic construction matters, as shown by changing it while holding the static action fixed. The quantitative mechanism is only partially resolved. Several registered spectrum checks failed. The post hoc contour checks corroborate many roots, but local contours around previously found roots do not establish completeness.

## GC34G

**Candidate screen.** The broader direction scan found that the per-cube least-squares fit grows in some directions and supports an extra slow wave. The per-vertex fit has an extra fast wave and growth at short wavelength. The blends still grow. All seven candidates kept the static kernel and passed the moving-source benchmark. Those checks alone cannot distinguish a stable tensor model.

**Candidate B.** Averaging the six tetrahedra’s strain rates before computing the kinetic energy is the best lead from the screen. Post hoc, no robust growing root or third low-frequency wave was found in the sampled cases. One polarization was about 5.8% fast in a sampled direction at kh = 0.05, against the 0.1% tolerance. The zero-frequency linear drift still needs classification as gauge or physical. The exploratory time runs are inconclusive. A smaller timestep reduced unseeded wave drift, but did not reliably reduce the apparent seeded growth. The code directly solves a timestep system built from a singular constrained formulation, then estimates the seeded components by subtracting trajectories and fitting their difference to modes. A constraint-reduced evolution with a direct tidal observable would give a better test.

## Priorities

1. Give B a bounded analytical follow-up. Identify the edge combination with no kinetic term, derive its constraint, and determine whether the linear drift changes tidal curvature. Use that reduction to count modes and evolve one Fourier component. Keep the wave-speed error as a separate gate.
2. If the drift is physical or the directional wave-speed error has no clear repair, stop tuning fits. Derive the next kinetic action around an explicit lattice identity preserving the freedom to change how time is sliced. Verify the identity and mode count before coupled-knot or longer runs.
3. Preserve the distinctions in the white paper. GC34D’s static result survives. GC34F confirms a defect in the time-dependent extension. GC34G gives requirements for a replacement. A stable replacement would still leave open how the original Model Gc rules could generate the shared geometry.

## Review scope and records

I reviewed the registered plans, result notes, source, and available evolution and check records. The full GC34G candidate spectrum was unavailable through the repository connector. I did not rerun the full spectra. The midpoint comparison was independently calculated from committed GC34F evolution data.

- [GC34F plan](../spec/EX-GC34F_plan.md), [results](../results/EX-GC34F/EX-GC34F_note.md), [evolution data](../results/EX-GC34F/evolve.json), and [source](../src/exgc34f.py)
- [GC34G plan](../spec/EX-GC34G_plan.md), [results](../results/EX-GC34G/EX-GC34G_note.md), [candidate B evolution](../results/EX-GC34G/evolve_posthoc_B.json), and [source](../src/exgc34g.py)
