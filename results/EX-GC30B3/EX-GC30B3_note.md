# EX-GC30B3 — Dynamical metric locking (EXPLORATORY)

**PI CLOSED 2026-11-22 (Richard: "I think those runs can be closed"; Claude review review/EX-GC30A-B3_claude_review.md).**

Interacting massless quantum fields can approach a common limiting speed dynamically, unlike the fixed tuning in GC30B/B2. Published one-loop Yukawa flow reproduces an infrared-attractive fixed line, but the coupling also weakens, so after an illustrative 40 decades in energy the initial ~20% speed discrepancy remains ~4% at alpha0=.03. Linearized multi-fermion QED shows a weak-sector bottleneck and exact disconnected neutral counterexample. Relative symmetric traceless tensor perturbations share the linearized attraction under stated operator-mixing assumptions, not a demonstrated robust full metric. The preregistered finite-scale gate fails. This excludes neither strongly interacting fixed points nor pregeometric theories, but no gravity sourcing should proceed on these benchmarks.

Design/run ChatGPT, PI approved 2026-11-22; exploratory, plan spec/EX-GC30B3_plan.md de961dca, LOG 79198042 before production. One numerical thread, 0.403 CPU s (2 core-hour cap). 20 full Yukawa trajectories, 35 linear gauge cases, 12 tensor cases. No later stage. All numerical checks pass; no validation rerun. Reference field theories only, not Model Gc.

## A. Yukawa flow: actual published one-loop equations

Use Anber & Donoghue (2011), https://arxiv.org/abs/1102.0789, Eq.18 massless scalar + Dirac fermion, rotationally invariant but independent limiting speeds b=cb,f=cf, coupling g. Their beta convention d/d ln mu. Implement beta_f = g²(f−b)/[6π² b(f+b)²], beta_b = g²(b²−f²)/[8π² b f³], beta_g = g³(3bf²+2b²f+b³+4f³)/[8π² b f³(f+b)²]. IR time t=ln(mu0/mu), integrate negatives. Ratio derivative (beta_b−(b/f)beta_f)/f agrees with Eq19 to 2.78e-17 and has sign(b/f−1) in UV beta convention; hence monotonic attraction in IR. At common c, beta_g=5g³/(16π²c³), independent velocity Jacobian eigenvalues 0 and 7g²/(24π²c³) (finite difference error 5.98e-11). Common velocity fixed line, not a unique velocity independent of initial parameters. With alpha=g²/(16π²c³) and eta=b/f−1 near unity, alpha'=-10 alpha², eta'=−(14/3)alpha eta. Integrating yields eta(t)=eta0(1+10 alpha0 t)^(-7/15); fixed alpha gives eta0 exp(−14 alpha0 t/3). This follows directly from Eq18; no fitted power. Compare paper discussion Sec7. The full integration initializes g from alpha0 referenced to c0=1 (not mean velocity); reported instantaneous alpha uses mean speed, so away from equality these differ. Alpha0=.1 is an optimistic marginal/strong loop stress not precision controlled; even .03 is only exploratory one-loop. No claims of physical scale calibration. 

| b0 (f0=1) | alpha0=.001 | .01 | .03 | .1 (stress) |
|---|---:|---:|---:|---:|
| .8, mismatch at 40 decades (%) | 15.575 | 6.984 | 4.331 | 2.510 |
| 1.2, mismatch at 40 decades (%) | 13.949 | 6.569 | 4.055 | 2.330 |

Mismatch defined as (max speed−min speed)/mean. At b0=.8, alpha0=.03, sequence at decades 0,1,5,10,20,30,40 is 22.222%,16.636%,10.312%,7.868%,5.875%,4.921%,4.331%; alpha drops to .00105716. At b0=1.2 alpha drops to .00103539. None of six 10–20%-scale initial large mismatch cases at alpha0<=.03 passes 1%; neither optimistic .1 passes. For small eta0=±.01 at alpha0=.03, final eta=−.00209123,+.00208876 vs linear ±.00208998 (relative ~0.06%). Matched-speed controls remain equal. Max integration refinement discrepancy 2.55e-9 (absolute state variable). Asymptotically attractive; finite scale criterion fails.

For initial |eta|=.1, linearized closed-form required decades D=[(.1/target)^(15/7)−1]/(10alpha0 ln 10). These are formal extrapolations of a massless one-loop benchmark, not physical predictions over such scales. At alpha0=.03:

| Target | Running alpha D | Hypothetical fixed alpha D |
|---|---:|---:|
| .01 | 199.702 | 7.143 |
| 1e-6 | 7.49806e10 | 35.714 |
| 1e-15 | 1.44765e30 | 100 |

1e-6 and 1e-15 are illustrative stringent targets, not specific observational bounds; no physical UV/IR hierarchy assigned. Constant-alpha control quantifies why a nonzero interacting fixed point could matter; not an actual four-dimensional Yukawa solution. The numerical conclusion follows because effective integrated interaction strength grows only logarithmically in t, not because the fields are uncoupled.

## B. QED network and weak-sector control

Linearize Eq33 of same paper near a common c with alpha=e²/(16π²c), charge factors q_i, one Dirac species per i. At leading order x_g'=alpha*(8/3)sum q_i²(x_i−x_g), x_i'=alpha*(16/3)q_i²(x_g−x_i), alpha'=-alpha²*(8/3)sum q_i² (standard one-loop QED running at common c). x=c_i/c0−1. Approximate isotropic near-common-speed model; not full nonlinear QED and no mixed Yukawa/gauge calculation. Contributions add at one loop; alpha running speed correction is higher order in velocity offsets. Matrix L has nonpositive eigenvalues, weighted invariant 2x_g+sum x_i; one zero mode when connected, two when one q=0. Exact linear matrix exponential time s=ln(1+b alpha0 t)/b, b=(8/3)sum q_i², verified by direct ODE to 6.74e-14. Initial x=(0,−.1,.1); max discrepancy relative to mean as above. Threshold diagnostic switches final edge off after 10 decades, retains its frozen matching coefficient for bookkeeping (not a light propagating massive particle prediction).

| q1,q2 | 40-decade mismatch alpha0=.03 (%) | threshold control (%) | constant alpha control (%) |
|---|---:|---:|---:|
| 1,1 | 1.271 | 2.978 | .00000796 |
| 1,.3 | 8.562 | 10.621 | 2.412 |
| 1,.1 | 12.541 | 12.870 | 10.907 |
| 1,.01 | 13.192 | 13.195 | 13.162 |
| 1,0 | 13.198 | 13.198 | 13.187 |

For q2=.1 the slow nonzero eigenvalue magnitude .07103 vs 5.333 for equal charges, about 75x weaker; q2=.01 magnitude .0007111; neutral it is zero. No guarantee a weakly interacting/neutral field is connected by some other interaction; electrically neutral matter can still lock via Yukawa or other coupling (not tested by disconnected control). An exactly decoupled kinetic sector is an independent geometric clock and disproves unconditional universality. The threshold is an illustrative piecewise leading-log diagnostic, not actual threshold matching or full matter EFT. At equal charges even this linear model falls just short of 1% under running; the finite-scale negative claim is much stronger for weak/disconnected sectors. No precision inference from 10% initial offsets in a linearized approximation. Source Roy, Juricic, Herbut (2015), https://arxiv.org/abs/1510.07650, shows interacting fermion/boson/gauge critical examples with common terminal velocity; see their distinctions 2 vs 3 spatial dimensions. Not a new derivation of those full beta functions.

## C. Symmetric traceless metric components (conditional)

At the Lorentz-invariant common-metric reference, constant symmetric traceless rank-2 perturbations can be represented by sources to scalar and fermion stress-tensor operators (modulo improvements/total derivatives, with consistent normalization). Lorentz covariance of the linearized mixing means the species-space map cannot depend on which component of that irreducible Lorentz tensor is chosen. Relative isotropic speed perturbation is one such component after removing trace/normalization and common coordinate change. Assuming this closed two-operator subspace and Lorentz-symmetric regularization, the IR species mixing at one loop is alpha*[[-2/3,2/3],[4,-4]], same as isotropic fractional speeds; invariant 6x_f+x_b, one common source zero mode, eigenvalue −14alpha/3. Applied independently to unit Frobenius-norm traceless diagonal isotropic, spatial x−y, xy shear and tx tilt sources with initial relative amplitude .01. All four give final amplitude .00208998 at alpha0=.03, after 40 decades. Common tensor shift remains invariant exactly. This is an analytic covariance inference and numerical check of the linear map, not an independent anisotropic loop calculation, not nonlinear tensor locking or birefringent gauge propagation. Constant source approximation excludes gradients, masses, other spin-2 operators, gauge constraints, and possible nonmetric Lorentz breaking operators. Without these assumptions the full tensor conclusion remains unresolved. Common conformal factor, massive clocks, and local collective deformation dynamics not determined. An attractive uniform cone is not an emergent energy-carrying massless graviton; Weinberg–Witten constraints remain relevant.

## Decision and next recommendation

The common limiting speed is dynamically attractive in connected tested models, a conditional positive C3 direction. The finite-scale universal gate fails (no large initial mismatch Y passes, weak/neutral G fails). Do not proceed to GC30C; do not spend additional precision on this slow-running reference. A nonzero IR interacting fixed point with large enough irrelevant relative-metric anomalous dimensions and a sufficiently broad scaling regime could offer a mechanism; require derivation rather than setting constant alpha by hand. Strong interactions may also destroy stable particle excitations or generate gaps, so this is not automatically an escape. A uniform Lorentz fixed point does not supply a dynamical shared metric, universal stress coupling, tensor pole, or G.

Recommend change of direction toward GC31A, a paper-only screen of pregeometric causal or quantum-correlation structures that might intrinsically define one operational causal order for distinct excitations; compare to interacting fixed-point reserve. Require independent retarded propagation and clocks, no inserted continuum geometry, causal order plus scale/volume, and tests against same static correlations but different dynamics. This is not a reason to presume entanglement transmits a force. Sakharov curvature stiffness remains a downstream question contingent on a shared dynamical geometry. PI decision required before further work.

## Reproducibility

Run `python src/exgc30b3.py`. Compact yukawa.json, gauge.json and tensor.json contain trajectories, analytic controls, scale estimates, fit/refinement and platform metadata. Same-author reproduction of published beta functions, not new loop calculation or independent external replication. Review/PI closure pending.
