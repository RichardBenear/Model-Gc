# EX-GC28B to EX-GC28G — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-22. Reviewed ChatGPT's EX-GC28B (651bb1b–566716b), 28C (cb12f8a–a44d6d2), 28D/E/F
(1c19bd1–3bad391) and 28G (7454253–0950b87), plus the response to my previous review (1ba416d). These are scalar QFT,
Skyrme and GR reference models, not Model Gc. Closure is the PI's.

## Verdict

All six are correctly computed and honestly reported, and every run reproduces. I found no mathematical error. The
headline conclusions are right, and most of them are known results, recovered here as calibrations:

- **28B:** the zero-point (Casimir) force between two localized defects in 3D falls as 1/d³. It is not Newtonian.
- **28C:** a symmetry-protected massless collective mode cannot carry a universal monopole source.
- **28D:** Sakharov's induced Einstein–Hilbert term, with the textbook coefficients.
- **28E/F:** Gauss's law and screening by a polarizable medium.
- **28G:** with minimal metric coupling, the whole dressed particle (binding included) gravitates through its ADM mass.

Two of the negative findings are partly built into the chosen couplings. 28G's "accounting solved" leaves out a
long-range force that its own model contains. Both points are below. Neither needs a rerun.

## What I checked

1. **Reproduction.** I ran every script on Linux (numpy 2.5.3, scipy 1.18.1): 28B (both stages), 28C, 28D/E/F and 28G.
   Compared with the committed records, 28C/E/F/G match exactly and 28B to 2e-15. 28D's only difference is a fitted
   coefficient that should be zero: 9.59e-11 vs 9.59e-11. Total CPU was about 2 minutes.
2. **28B against theory.** The weak-coupling 3D result is E = −Q₁Q₂/(64π³d³) = −a₁a₂/(4πd³), with a = Q/4π. I
   derived it from the leading term of the TGTG formula, (1/2π)∫dξ Tr(T₁GT₂G) with G = e^{−ξr}/4πr. The records agree
   to 0.02–0.6%. The 1D opaque limit −π/(24d) is the standard result. This is the outcome I suggested registering in
   advance, and the plan did state it before the run.
3. **28D against analytic values.** The values are a₁ = (1/6 − ξ)R per scalar and −R/12 per Dirac component (−R/3
   total). At Λ = 10, M = 1: ∫e^{−t}t^{−2}dt = 94.967, and /16π² = 0.6014 gives 0.1002 per minimal scalar and 0.2005 per
   Dirac field. That matches the table exactly. The signs (minimal scalars and fermions positive, conformal zero,
   ξ > 1/6 negative) are the known Sakharov results.
4. **28G field equations.** I checked m′ = dE/dr in the normalized units. δ′ = G(2Kf′² + r²φ′²)/(3πr) is 4πGr(ρ+p_r)/N
   after the 1/(12π²r²) normalization. I also checked the radial and tangential stresses, including the split of the
   quartic Skyrme term into (r,θ), (r,φ) and (θ,φ) pieces. Komar = ADM is then the Tolman/virial theorem, so agreement
   to 2e-7 (massive) checks the solver, as the note says.
5. **Cross-consistency.** 28E's undressed energy (1.26251) and 28G's flat neutral case (1.262506) equal 28A's radial
   classical energy with the 27B potential. 28E's translational mass formula (2E₂ + 4E₄ + 2E_field)/3 is the right
   von Laue combination: potential-type terms weigh 0 in it and 3 in the virial.
6. **Protocol.** Each study has a plan and a LOG row committed before its results commit. Each LOG row quotes a
   separate go-ahead from you. Each was well under its 2 core-h cap: the whole series used about 3 CPU-minutes, and
   everything is labelled EXPLORATORY. Two deviations were registered and disclosed, with the original rows kept: the
   28B infrared quadrature fix and the 28G restart. The 28A response accepts my earlier points and adds the
   Born-derivative table I asked for. Its finding that the subtracted term keeps growing with the cutoff, and that a
   comoving cutoff changes it, is the honest answer. One caveat: plan-to-result gaps are 1–3 minutes, so registration
   is formally before the run but not a separate step anyone could review.

## Issues

**1. 28C Arm A: the "universality fails" result is built into the coupling (substantive).** The collective scalar
couples as M² ± gs, so it is odd: by construction it couples to the difference between the two species and never to
their sum. A common positive shift (λ₊ = λ₋) then gives J = 0 by symmetry. That is the coupling choice, not a finding
about vacuum response. The natural even, dilaton-like alternative, M²(1 + gs) for both species, couples to the sum
(the trace of the mass matrix). It would give a source of the same sign for every positive-energy defect, and it was
not tested. So the note's "natural defect sources vanish or reverse sign" holds for this coupling only.

**2. 28C Arm B is a theorem, and that is its real lesson (framing).** The protected massless angle is a Goldstone
mode. The shift symmetry that keeps it massless also forces it to couple only through derivatives, so a static defect
cannot give it a monopole. That is what the pinned and free-rotation tests show. The general statement is that
symmetry-protected masslessness and monopole sourcing work against each other for a scalar. That is a known obstacle to
scalar gravity, and the reason gravity's own protection comes from gauge (diffeomorphism) symmetry. Saying this
directly would make 28C's conclusion general instead of model-specific.

**3. 28E/G: the secondary scalar is a second long-range force, and it is not universal (substantive).** 28G shows the
metric sees total energy whatever the secondary charge Q is. But the secondary field is still massless, so two dressed
particles also exchange a 1/r scalar force ∝ Q₁Q₂, attractive for like charges. Within this cohort Q/M changes about
5× with α at almost fixed mass, and flips sign with q. In an equivalence-principle test that force counts, so the
combined attraction would depend on composition unless Q/M is tiny or constant. "Separating the signed source from the
metric source solves the accounting problem" is right for the metric. The full model, though, adds a fifth force that
universality forbids. I'd add one sentence to the 28G note.

**4. 28G: the universality comes from the assumption (framing, already stated).** The note says correctly that minimal
metric coupling is an input. Read together, 27A/B, 28D and 28G bring the programme back to where EX-GC25/26 left it:
coupling to the metric gives universal sourcing, and the open problem is why matter couples to a common metric. 28G
does improve on the 26 series in one way. It uses full GR (tensor), not Nordström scalar gravity, so the light-bending
objection to the 26 setup does not apply to it.

**5. Minor.**
- 28C: the J/E ratios compare a log-divergent source with a quadratically divergent bare energy. The note already flags
  this as not physical, and its fall from .504 to .0138 between Λ = 2 and 20 shows why.
- 28E/F/G use ρ = (1 − cos f)² as the secondary source. That is a free choice, so the specific α/q numbers carry no
  meaning outside this model.
- 28G: the two failed tighter-tolerance massless runs are disclosed. The wider-domain controls converge, so I see no
  concern.
- 28D: "Ward check" is a check of the linearized Einstein operator's algebra, as the note says, not of the determinant.

## Bottom line for the PI

- The numbers are right and reproduce.
- What the series establishes is mostly negative and mostly known:
  - Zero-point Casimir forces are too short-range (1/d³).
  - Vacuum fluctuations can induce a gravitational stiffness (Sakharov), but the cutoff fixes its size.
  - A protected massless scalar cannot be sourced like gravity.
  - Universal sourcing appears only once metric coupling is assumed.
- I'd count 28B and 28D as clean closures.
- 28C's negative conclusion is weaker than written for Arm A (issue 1).
- 28E–G should mention the secondary scalar's own non-universal long-range force (issue 3).
- The question the series leaves open (where common metric coupling comes from) is the same one EX-GC25/26 left open.
  Further loop calculations in assumed backgrounds are unlikely to answer it.
