# EX-GC28C — loop-induced collective response (EXPLORATORY)

**Summary:** In a Lorentz-covariant scalar-field toy model, integrating out microscopic vacuum fluctuations generates a positive kinetic stiffness for a collective variable that has no bare kinetic term. A fine-tuned ordinary scalar and an exactly rotation-symmetric angle both support a 3D 1/r static response to an externally supplied monopole; the latter needs no mass tuning at the one-loop quadratic level. However, naturally generated weak static defect sources either vanish or change sign with internal response, even for equal positive Born energies. Thus induced propagation succeeds, universal gravitational sourcing does not. This is a useful separation of the problems, not an emergent gravity claim or a test of the numerical strength of real gravity.

Design, implementation, run: ChatGPT. Reviewed (review/EX-GC28B-G_claude_review.md); PI CLOSED 2026-11-22 (Richard: "go ahead and close 28 and 29A and we will defer 29B"). Plan `spec/EX-GC28C_plan.md` and LOG committed first. Source `src/exgc28c.py`, run `python src/exgc28c.py`. 24 model/detuning combinations, 144 static radii (plus loop-source probes in scalar arm), 28 defect source rows, 144 kernel samples, 8 validation rows. 2.02 CPU-s, one thread, within 2 core-hour cap. No deviations.

## Motivation from PI's Puthoff comment

H.E. Puthoff, *Gravity as a zero-point-fluctuation force*, PRA 39 2333 (1989), https://doi.org/10.1103/PhysRevA.39.2333, proposed EM ZPF-induced particle Zitterbewegung and gravitational interaction. S. Carlip, https://doi.org/10.1103/PhysRevA.47.3452, identified a calculation error and reported a negligible nongravitational force; Puthoff's reply https://doi.org/10.1103/PhysRevA.47.3454 argued physical cutoffs restore inverse square. The scalar loop test here is inspired by a secondary response rather than direct vacuum force but does not reproduce that disputed electromagnetic calculation. Sakharov-style induced effective action is a different mechanism (see M. Visser arXiv:gr-qc/0204062): a kinetic term for a collective field induced by quantum fluctuations, not ordinary Casimir attraction. We do not assume a gravity mediator is already dynamical.

## Model and quantitative result

Euclidean 4D effective action is S_local[s]+(1/2)Tr log D(s); microscopic real scalars have positive squared masses a,b. Auxiliary s has no bare gradient term. At one loop, with Lorentz-covariant proper-time cutoff Lambda, bubble B_ab(z) = (16pi²)^-1 int_0^1 E1((xa+(1-x)b+x(1-x)z)/Lambda²) dx, z=p_E². Induced kernel has positive slope Z because B decreases with z. Bubble identity with tadpole, finite-difference slope, rotation invariance verified to 1.06e-11, 1.73e-13, and machine precision respectively. The angle seagull and bubble cancel at zero momentum. No-loop local kernel is momentum independent (contact response smeared only by source size), or zero/undefined response for a completely flat auxiliary angle. Thus the nonlocal kinetic stiffness is induced, not inserted. Lorentzian pole continuation has positive residue and omega²=k²+m_pole² below microscopic two-particle threshold. Lorentz invariance of the underlying continuum model and cutoff is stipulated, not derived.

**Arm A: ordinary auxiliary scalar.** Microscopic squared masses 1±s, local potential u s²/2+s⁴/4. K(z)=u-B_11(z). For u=0, K(0)<0: the reference vacuum is unstable, not a massless gravitational field. At ucrit=B(0) the quadratic kernel is gapless and local quartic is positive (~.99895; local stability only, not an asserted global UV completion). Detune by delta=Z mu², and the pole mass is close to mu. This is critical tuning, not a protected long-range mode.

**Arm B: symmetry-protected induced orientation.** R(s) diag(1,4) R^T(s). Constant angle just rotates the microscopic species, hence determinant independent of s. K(z)=(a-b)²[B_ab(0)-B_ab(z)]. Zero mass does not require adjusting a local mass counterterm; exact global rotation symmetry forbids angle potential. Degenerate a=b control Z=0 (no physical relative orientation). In nondegenerate model the one-loop effective kinetic term is positive. The result is conditional on the ordered background and quadratic approximation; quantum fluctuations of s, phase disorder, nonlinear defects and full dynamical completion are not tested. This is a Goldstone-like coordinate, not a spin-2 metric.

| cutoff Lambda | scalar Z | scalar ucrit | angle Z |
|---:|---:|---:|---:|
| 2 | .000821969 | .00661300 | .00235352 |
| 5 | .00101405 | .0169793 | .00376474 |
| 10 | .00104493 | .0255705 | .00403254 |
| 20 | .00105279 | .0343020 | .00410268 |

Both stiffnesses approach finite positive values; scalar mass cancellation is cutoff dependent. The angular symmetry protects the constant mode for every cutoff used. With unit Gaussian probe, Fourier sine integration of the full kernel at Lambda=10:

| r | scalar tuned | scalar mu=.05 | angle protected | angle mu=.05 |
|---:|---:|---:|---:|---:|
| 4 |19.0420|15.6028|4.93347|4.04337|
| 16 |4.75975|2.14050|1.23336|.554779|
| 64 |1.18994|.0485602|.308341|.0125838|
| 128 |.594969|.000990111|.154171|.000256518|

At r>=16 the tuned/protected response matches 1/(4pi Z r) within numerical precision. Detuned response matches the pole Yukawa form with positive residue and a finite Gaussian form factor; at Lambda10 mu=.05 masses .0499937/.0499972. For mu=.2 masses .199597/.199820. Suppression is not a change in numerical resolution. The unit probe is imposed and is **not** a natural gravitational monopole. An induced loop source with Gaussian profile has the same asymptotic tail when nonzero (normalized to J(0)=1 in records; physical amplitude is J(0) from source table). Attraction for two static same-sign sources follows at quadratic level from -J1 K^-1 J2, repulsion for opposite signs.

## Source test: the blocking issue

Arm A weak localized static positive squared-mass perturbations lambda_+,lambda_- yield cutoff-regulated first-order vacuum energy E=(lambda_++lambda_-) I(1)/2 and monopole J=(lambda_+-lambda_-) B_11(0)/2. A source from the loop (secondary field), not an imposed bare term. At Lambda=10:

| perturbation (lambda+,lambda-) | E | J | J/E |
|---|---:|---:|---:|
| (1,0) | .300693 | +.0127852 | +.0425193 |
| (0,1) | .300693 | -.0127852 | -.0425193 |
| (1,1) | .601386 | 0 | 0 |
| (2,0) | .601386 | +.0255705 | +.0425193 |

First two are related by exact interchange symmetry; same energy, opposite source. Next two have equal Born energy and distinct sources. Common species shift cannot generate an odd scalar source by symmetry, even beyond first order at s=0. J/E for one species changes .50426→.0137626 between Lambda2→20; bare cutoff energy is not physical renormalized total mass. Even a positive common energy does not enforce universal sourcing. Defects are external perturbations with infinitesimal strength parameter; no inertial masses or nonlinear backreaction calculated.

Arm B equal diagonal perturbations preserve vanishing angle torque. Adding opposite pinned off-diagonal elements +/-.5 gives same positive Born energies .570948 but J=-/+ .0304380 at Lambda10. This explicitly fixes an anisotropy relative to the vacuum, not an invariant mass. If defect can rotate with the background, constant angle symmetry is restored and uniform derivative of energy is zero. A relaxed orientation cannot sustain a net static torque without external pinning; derivative/multipole sources are not universal energy monopoles. The exact protected pole is therefore insufficient.

## Decision and next step

Positive: vacuum fluctuations can generate stable-sign kinetic stiffness, with a massless pole protected by symmetry rather than tuning in an ordered one-loop model. A genuine secondary long-range 1/r response is possible, distinct from 28B's two-propagator Casimir loop. Negative: none of the natural weak source tests yield universal energy sourcing. The coefficient is in toy units, depends on chosen coupling, symmetry and regulator; no prediction of Newton's G (or Puthoff's magnitude claim). The microscopic variables and Lorentz symmetry were chosen; not full induced GR, no metric or spin-2 gauge dynamics, no dynamical topology, no guarantee of nonperturbative order. No inference about actual charged-particle Zitterbewegung.

Recommend EX-GC28D as an analytically controlled Sakharov metric-action screen: compare the stress-tensor/metric effective action from integrating out scalar and fermion vacuum fields; check diffeomorphism Ward identities, induced Einstein-Hilbert coefficient and sign, the cosmological term, cutoff dependence and whether universal energy sourcing follows from the assumed metric coupling (not call that emergence). Include nonmetric control and check what protects the gapless pole. A separate Puthoff-inspired oscillator/EM response experiment could test secondary-field power, phase averaging, retardation and neutrality without assuming gravity, but must confront Carlip's critique and avoid fitting cutoffs to G. The main bottleneck is now a protected collective mode with universal sourcing, not attraction alone.
