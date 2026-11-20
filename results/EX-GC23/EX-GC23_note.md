# EX-GC23 — advected velocity and actual particle coasting (EXPLORATORY)

PI Richard Benear; design, code, runs and analysis: ChatGPT; motivation and review of EX-GC22: Claude. Our construction, not LaViolette's. Based on f4b7198 (EX-GC22 review merged). Pre-run remote plan commit fc1f500; `spec/EX-GC23_plan.md` and `review/EX-GC23_response.md`. No confirmatory claim or PI closure.

**Summary.** The material velocity derivative repairs most of the earlier artificial loss of motion. At h=.25 the W centroid now closely follows the flow, and the strong particle travels 0.96345 of the ideal 1.0 after a 20-unit kick while retaining 97.06% of its internal speed (uncorrected 73.91%). The weak particle retains 99.81%. The momentum ledger identifies the dynamic cancellation of the consumption-reweighting term. Both types still fall at almost the same actual acceleration under equal applied flux. The remaining slow decay is not zero, the momentum density is assigned rather than derived, and the distinct EX-GC17 head-on absorption drag is disabled. No emergent Newtonian inertia or Lorentz invariance is established.

## 1. Model and controls

The EX-GC22 Gc1 chemistry and wave medium are unchanged (no knot, D=80, c=2, tauJ=20, h=.25, L=6 except wall control). The chemical deviations from ambient are advected by V; G remains the independent medium. `R=wgt*(2J/G)/20` with `wgt=s/(s+0.001 max(s))`, `s=KGW`. The original cohesion C is the local s-weighted Gaussian average of V minus V (sigma=1.5). The transport switch changes only the equation for V:

```
off: dV/dt = R+C
on:  dV/dt = R+C - (V dot grad)V
```

Matched independent 30-unit Eulerian rests, stationary G from Poisson iteration, then correctly initialized face flux J=-80 grad G. A uniform +.05 boost is applied to the full V (no initial chemical shift). Flux trials add uniform Jz=.02 with E=-.02/80 to G and J from t=0. No head-on drag. dt=.000375, velocity step .0075, same reflecting chemistry and Dirichlet medium ghosts as EX-GC22. Four workers at most, numerical libraries single threaded. All rest preparations use original EX-GC22 chemistry stencil; vectorized algebraically equivalent version used for production, validated against both original and NumPy.

## 2. Matched kick results (20 units, L=6)

| K | Transport | W shift | V(20)/V(0) | Late W speed | Late gamma from V | Max fractional M change |
|---:|---|---:|---:|---:|---:|---:|
| .2 | off | .847354 | .739108 | .0389112 | .0160150 | .000479568 |
| .2 | on | .963450 | .970598 | .0477535 | .00219527 | .000404482 |
| .01 | off | .967514 | .982927 | .0480972 | .000913137 | .000619011 |
| .01 | on | .974762 | .998139 | .0486594 | .000140931 | .000609231 |

Nominal ballistic displacement is .999994 (rounds to 1); late means the fit interval t=10–20. Strong speed loss is reduced by 88.7%; late decay rate by 86.3%. Actual late W speed / mean flow is .974 to .978 in all four runs (coarser EX-GC22 range .38–.58). This reduces the lattice-slip confound, but h=.25 is one resolution and does not establish a continuum limit. The remaining W speed deficit persists even at weak K. All four cores survive, no M changes exceed .07%, W peak shifts <1%. s centroid agrees with W centroid to <.001 model unit. No indefinitely maintained coasting shown.

## 3. Exact discrete momentum accounting

Define P=∫sV (omit common assigned tau=20). Between velocity sweeps, V is fixed; changes in P are exactly changes of s at fixed V (density changes due to all chemistry, medium transport and chemistry advection, not only local creation). In a velocity sweep s is fixed, and increments from recoil, cohesion and material transport are accumulated separately. At every sample the sum is checked against P(t)-P(0). No claim that P is conserved; the purpose is to identify contributions.

| K=.2, kick | Recoil | Cohesion | Transport | Changed s | Net delta Pz |
|---|---:|---:|---:|---:|---:|
| off | -.00924850 | -.00629284 | 0 | -.0646993 | -.0802406 |
| on | -.0103206 | +.0000587311 | +.0750751 | -.0737475 | -.00893430 |

Transport balances the changing density contribution to about 1.8%, dynamically confirming EX-GC22. Cohesion's integrated force also becomes much smaller when the pattern and internal V travel together. The dominant remaining negative ledger component is recoil, but not necessarily the steady wake's entire magnitude: finite-box transients, evolving source and medium relaxation at tauJ=20 are present. The EX-GC22 prescribed Gaussian gamma (~.000235 for K=.2 D=80) cannot be directly subtracted from this late rate. The residual closure errors for the four kick runs are below 5e-13 absolute. This does not establish physical total momentum; in particular s is created/destroyed by the chemistry and the dynamics of the recoil partner has not been accounted for.

## 4. Imposed flux response (transport on)

| K | Settled M | Flow accel [0,10] | W accel [0,10] | W shift at t=20 |
|---:|---:|---:|---:|---:|
| .2 | 6.15961 | .00198472 | .00192033 | .384977 |
| .4 | 8.71817 | .00198032 | .00191351 | .383067 |

Nominal reference 2J/(G tau) ~.002 (G about 1, weighting and internal self-push affect it). Flow accelerations differ by .22% and actual centroid accelerations by .36% relative to their mean; actual early acceleration is ~96% of nominal. Both cores survive (<.05% M variation). This is near-type-independence only under the assigned mass-per-absorption rule and equal imposed flux. It is not an independent test of the equivalence principle, nor a two-body interaction or orbit.

## 5. Null and conservative alternative

The K=.2 transported null stays stationary throughout 20 units: |W centroid| at end <3e-17, |Vcore,z| <4e-18. M changes -.01595%; inward radial V survives. Turning transport on does not cause an otherwise stationary particle to self-accelerate.

The conservative alternative uses mass-advection flux `dm=-div(sV)` and momentum-advection flux `dp=-div(sVV)`, same sealed-face upwind donor, and velocity increment `(dp - V dm)/s`; the actual chemical s still evolves from the independent equations. Thus it does **not** give an independent conserved inertial density. The symmetric cohesion exchange `C=[G_sigma*(sV)-V(G_sigma*s)]/max(G_sigma*s)` has zero integrated s-weighted net force by symmetry of the convolution; unlike the original per-site normalization it is weaker in dilute regions. No extra tuning. The alternative is compared from the original settled start, not a newly equilibrated conservative-cohesion rest state.

At t=20: W displacement .961259, V retention .965098, late W speed .0475777, late gamma .00246326. These are close to the material-derivative result, so the strong coasting improvement is not reliant on nonconservative Gaussian cohesion. Consumption falls 4.87%, W peak changes <1%, survival threshold passes. Momentum ledger (z): recoil -.00953115, cohesion +3.36e-19, velocity transport +.0980530, changing s -.113746, net -.0252243. The conservation of cohesion is verified, but loss of density still means the assigned particle momentum is not conserved. The altered cohesion also compresses the core; a later conservative closure would need its own equilibrated rest and an independent continuity/source law if physical inertia is sought.

## 6. Wall control

The ±9 box (73³ instead of 49³) has its own 30-unit settled rest; M=6.16242 vs 6.15961 in ±6. Budget adjustment before the wall kick: T=10 instead of 20, enough for about one round-trip wave flight at c=2 (9 units), not enough to distinguish later reflections from long-time relaxation. Matched time samples:

| t | ±6 W shift | ±9 W shift | ±6 Vcore | ±9 Vcore |
|---:|---:|---:|---:|---:|
| 4.99875 | .243978 | .243987 | .0498379 | .0498396 |
| 9.99750 | .485511 | .486092 | .0495898 | .0497512 |

At t≈10 W shifts differ by .120%, Vcore by .325%. The larger box lowers the early gamma from .000832255 to .000221213; residual loss is sensitive to geometry/relaxation even though the absolute difference is small. G depletion differs by ~9.8% at t≈10. The main transport improvement is robust to this modest box-size change over the measured interval; the small remaining wake/recoil decay cannot be claimed independent of boundaries. At terminal t=10.000125, larger-box M changes -.0238% and V retention .995076.

## 7. Interpretation and recommended next step

The material derivative is a kinematic consistency requirement for the advected pattern; it is not itself a microscopic derivation of inertia. The paired trial shows that the main previous drag-like signal was an Eulerian transport defect, and a fine lattice makes W track V. The surviving smaller recoil signal is measurable, not to be rounded away. The bare absorption drag, a separate medium-rest-frame effect, is still deliberately absent. The G telegraph equation with finite relaxation time and the chemistry have not been shown to be Lorentz invariant; the model's explicit flow degree of freedom and effective mass proportional to s were introduced by construction.

**Candidate EX-GC24 (not approved):** confront the explicit EX-GC17 head-on absorption drag with the now corrected flow and fine-grid centroid. First reinstate the original drag term in matched strong/weak kicks and flux trials; check decay relative to the expected assigned relaxation time, distinguish it from the residual self-wake recoil, and quantify the corresponding centroid response. Then test a separately specified momentum-exchange law that accounts for recoil of both pattern and medium, using a conserved or explicit sourced momentum ledger. Uniform relative motion and a uniform external gradient should be distinct controls. Do not remove a real drag term merely to fit Newtonian motion, or infer an equivalence principle from the assigned absorption-proportional inertia. A two-body orbit test is premature until this separate obstacle is addressed. A h=.125 subset would quantify residual lattice slip, and longer wave runs would resolve finite relaxation transients before extrapolating to steady coasting.

No result here is a direct prediction for a Standard Model particle, quantum gravity, or a physical measurement; this remains an exploratory model result. No claim of indefinite coasting, physical energy, total momentum conservation, or inertia derivation.

## 8. Reproducibility and compute ledger

Nine production runs (four matched kicks, two flux, one null, one conservative, one larger box) and four independent rest preparations completed. All production cases have status complete. Sum of measured CPU time 5815.96 s = 1.61555 core-hours, plus <.01 core-hour for brief validation, benchmarks and reporting, within the 2 core-hour budget. Maximum 4 concurrent workers; deterministic (no random sampling in production). Platform, package versions and individual wall/CPU times in the records. Main protocol deviation: larger-box T=10 rather than 20 (prospectively logged for the wall arm); accelerated stencil before kicks and added progress saves, also logged. No excluded outcomes.

- `python src/exgc23_transport.py --phase validate` (six checks; chemistry, wave, affine transport and conservative ledger operators). `python tests/test_regression.py` passes. The fused wave step agrees with inherited NumPy exactly on the validation fixture; original chemistry vs NumPy <3.5e-13; fast path vs original <6e-14; other tests <4e-16. Max momentum ledger residual across the nine runs <5e-13. Linux C++ compiler required for this implementation.
- `python src/exgc23_transport.py --phase settle --workers 4` regenerates four checkpoints (original EX-GC22 stencil); `python src/exgc23_transport.py --phase run --workers 4` runs the listed finite cohort (production optimized stencil). For matching the exact execution schedule use case names via `--cases` and see plan. The initial arrays live in `local_raw/EX-GC23/settle_*.npz` (four files, hashes embedded in their JSON records); separately provided `EX-GC23_settled_fields.zip` contains those exact archives, verified against all four SHA-256 hashes. They are not in git. Extract at repo root to reproduce production without settling.
- `python src/exgc23_summary.py` recomputes all key fits, ratios, ledger residuals and compute totals from the compact JSON without solver import. `python src/exgc23_figure.py` generates `paper/figures/exgc23_transport.png`. Figure and summary regenerated after the last control. The review disposition is in `review/EX-GC23_response.md`.
- Branch `ex-gc23-chatgpt` isolated from main; PI assigns review/merge and any follow-up stage. No external independent replication.
