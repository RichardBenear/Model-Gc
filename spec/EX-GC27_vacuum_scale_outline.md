# EX-GC27 — vacuum-scale sensitivity (draft outline, not registered)

PI direction, 2026-11-21: prioritize Candidates 2 and 4 in [Candidate Models](../paper/Candidate_Models_next_phase.md), beginning with the inexpensive vacuum-scale sensitivity test. Claude is separately responsible for EX-GC26E3. No production authorization implied by this outline.

## Question

Can the source of gravity be identified with the response of localized particle energy to a property of the vacuum? Can that response be universal without assigning each particle its own gravitational coupling? Start with existing conservative Skyrme reference states (not Model Gc), rather than building a new particle solver before testing the mechanism.

Use c=1, a dimensionless vacuum coordinate eta = ln(chi/chi0), and define, for fixed topological sector and branch,

$$s_i=\left.\frac{d E_i}{d\eta_\infty}\right|_0,\qquad k_i=\frac{s_i}{E_i(0)}.$$

E is excess energy above the corresponding homogeneous vacuum. A positive common k is a candidate universal scalar-source coupling, not evidence by itself of gravity. For a later canonically normalized scalar with energy f_chi^2 (grad eta)^2/2 and a local matter action, a static weak source s produces eta ~ -s/(4 pi f_chi^2 r); the reciprocal force is -s_test grad eta. These statements require a massless scalar, a stable equilibrium, and consistent normalization. Neither masslessness nor the stiffness f_chi is derived here.

## A crucial analytic control (do before any production)

Ordinary massless Skyrme energy E = A E2 + B E4. Under a spatial dilation by lambda (size increases by lambda), E2 -> lambda E2 and E4 -> E4/lambda. For A = exp(a eta), B = exp(b eta), an exactly relaxed continuum branch has

$$E_i(\eta)=E_i(0)e^{(a+b)\eta/2},\qquad R_i(\eta)=R_i(0)e^{(b-a)\eta/2}.$$

This holds for every static branch if there is no other length scale. The envelope derivative is s = a E2 + b E4 at eta=0; Derrick balance E2=E4 makes k=(a+b)/2. Thus the simple quadratic-only coupling (a,b)=(2,0) is exactly the Nordström static result already measured in EX-GC26, not an independent discovery. A simultaneous common rescaling (2,2) is also imposed universality. The size-only control (2,-2) predicts zero source despite a changing core radius; it separates vacuum distortion from gravity.

## Phase A: inexpensive static screen

Use the hash-verified B=1 and B=2 rest states from 26E1/E2 at h=0.2 and 0.15, common L=7, matched discretization. Include vacuum and, if cheap, the parity-reversed B=-1 state; topology and branch identity fixed before revealing sensitivity. No force, no time evolution, no need to wait for 26E3.

1. Compute baseline E2, E4, E and core radius, geometric degree, stationarity residual and boundary sensitivity. For eta ±0.02 and ±0.05, evaluate energy at fixed shape, then at the analytically rescaled shape. Compare analytic envelope prediction to centered finite differences.
2. A limited set of true relaxed continuations (coarse first) tests whether shape response and numerical errors follow the dilation prediction. Reuse baseline fields; warm-start at scaled radius; retain all failed relaxations. Report residual/energy change, not merely iteration count.
3. Run couplings (a,b) = (2,0), (0,2), (2,-2), (2,2). All four are controls, not new candidate physics. At h=0.2/0.15 the inherited Derrick imbalance creates different finite-grid sensitivities; report it directly rather than renormalizing each species.
4. Measure the binding sensitivity s_bind = 2 s_B1 - s_B2 and compare with k E_bind. If E_bind is small, show absolute error as well as its relative value. Parity should leave energy and scalar source unchanged.

## Phase B: break the trivial scaling identity (the real candidate screen)

Introduce an independently motivated second vacuum scale. Minimal concrete test: a nonnegative potential V0 = integral (1-phi0) d^3x with positive coefficient C0, so E = A E2 + B E4 + C V0. Set C0 once using a preregistered modest mass scale (e.g. m R_B1 ~ 0.3; define normalization precisely in plan). Use the same C0 for all species; relax new baselines before derivative measurement. For C=C0 exp(c eta), the envelope source is

$$s=a E2+b E4+c E_V.$$

Derrick balance becomes E2-E4+3 E_V=0. With (a,b,c)=(2,0,4), the potential's coefficient scales consistently with the original dilation, and k=1 still follows by construction. With (2,0,0), k = (E-4 E_V)/E and generic species differences can appear. With (0,0,1), k=E_V/E is a nonuniversal sensitivity control. The point is to identify which independent terms/fractions break universality, not retune exponents until k matches.

Compute all sensitivities analytically from each stationary baseline, with centered finite-difference checks at eta ±0.02; only relax the subset required to check the prediction. A configuration whose k happens to agree with B1 at one resolution is not a discovery; compare independent energy partitions, sizes and additional branches. Do not change couplings after inspecting species outcomes. If the available states cannot discriminate because their partitions are too similar, record this as unresolved and add a structurally distinct particle in a separately approved phase.

## Interpretation and budget

- A: recover known scale identity and numerically characterize lattice error. Success is a calibration, not discovery.
- B: if distinct species sensitivities diverge when the second scale is varied independently, this is positive information: scale symmetry or a stronger microscopic principle is necessary for universal sourcing. If equality survives despite a nontrivial independent perturbation, check whether another exact symmetry or accidental equal energy fractions explains it before extending.
- Do not interpret equal k as evidence for a massless mode, Newton's constant, actual gravity, or quantum effects. A viable candidate must later source a dynamical vacuum response that is causal, long-range, conserves energy/momentum and coasts without intrinsic drag. An ordinary massive scale field gives Yukawa screening.
- Proposed initial cap 2 core-hours, with an explicit LOG row and frozen machine-readable plan before any production. Most analytic/fixed-shape measurements are cheap; restrict actual re-relaxations to a budgeted coarse subset first, then fine checks. If Phase B baselines need more time, report checkpoint and request approval rather than exceed cap. Keep all records and code in repo.
- This is exploratory, not a precision gate. Report k and its measured differences, finite-difference consistency, h/box dependence, rest-energy/size changes, degree and relaxation residual. No orbital work.

## Next step if warranted: connect Candidate 2 to Candidate 4

Sakharov's idea motivates asking whether the effective vacuum stiffness and coupling arise from quantum fluctuations, not merely assigning them. For a stable relaxed configuration, form the tangent-space quadratic fluctuation operator and kinetic metric, remove translational/internal zero modes explicitly, check negative eigenvalues, and compare generalized eigenfrequencies against the vacuum. Start with a finite, small lattice/1D control where independent diagonalization is possible.

Measure vacuum-subtracted spectral energy and the change of its derivative with eta. For two defects, evaluate the connected interaction E12-E1-E2+E0 at several separations; distinguish classical energy from the quantum correction. Vary cutoff and box, apply consistent counterterms and avoid interpreting a cutoff-dependent zero-point sum as physical source. A generic Casimir-like interaction is not Newtonian; measure sign, range and energy scaling rather than assume them.

Later, the Gaussian covariance of these same modes enables the entanglement/correlation branch, without equating correlation with force. The central open question is whether coarse-grained vacuum response produces a gapless universal coupling (eventually tensor geometry), not whether we can insert an energy-coupled scalar a second time.

**Decision requested:** approve Phase A and a bounded Phase B screen, or specify a different initial budget. No computation run and no production plan frozen yet.
