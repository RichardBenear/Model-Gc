# EX-GC28G — metric gravity of a dressed particle (EXPLORATORY)

**Summary:** Minimal metric coupling of the radial Skyrmion and secondary scalar yields a self-consistent gravitational field sourced by the complete dressed system, even as the secondary scalar monopole is reduced by screening, changes sign, or vanishes. ADM mass, the integrated matter/field mass constraint, and asymptotic lapse agree; an independent stress/Komar integral agrees up to the known finite-domain wall stress. Gravitational self-binding reduces ADM energy. This calibrates assumed generally covariant coupling and does not derive metric geometry, universality or G from the ZPF.

Design/implementation/run: ChatGPT. PI closure and Claude review pending. Plan `spec/EX-GC28G_plan.md`, registered pre-run; refinement failure and continuation remedy `spec/EX-GC28G_refinement_check.md`. Source `python src/exgc28g.py`. 80 baseline cases and 20 domain/tolerance controls, 98 converged, two tighter massless-core R24 runs (q=.25, alpha=0 or 5, G=.03, tolerance 5e-8) exceeded 30,000 mesh nodes, preserved as failures; standard-tolerance and R48 counterparts converge. Initial run stopped preparing a wide massless-core control after all 80 baseline cases had converged; final records were not written then, so repeated deterministic baseline, used converged metric seeds for refinements, checkpointed each record. Final run 53.33 CPU-s; initial run estimated under 2 CPU-min from execution duration, total under 3 CPU-min, well below 2 core-h, one thread. No change to cohort or criteria.

## Method

Same dressed hedgehog as GC28E, C=0 or .178396442851, q=0, +/- .25 and .5, polarization alpha=0,1,5,-.25 as cohort. Spherical static metric ds²=-exp(2delta) N dt²+dr²/N+r²dOmega², N=1-2Gm/r. Covariant Skyrme and secondary scalar action gives mass constraint m'=(N K f'^2+N r² phi'^2/2+U)/(3pi), delta'=G(2K f'^2+r² phi'^2)/(3pi r). Solve six coupled first-order BVP variables, regular at r=1e-4 and f(R)=0; exact massless scalar exterior integrated in x=1/r to zero, enforcing phi(infinity)=0 and normalizing delta(infinity)=0. This includes exterior scalar self-energy and metric backreaction rather than assuming a Schwarzschild exterior. The residual finite-radius f'(R) represents a truncated Skyrme tail and wall stress; the exterior is exact for the scalar-only continuation, not exact for the untruncated soliton. No horizon found (min N=.953873), no full coupled-mode stability analysis.

For local normalized radial densities A=N r² f'^2, B=sin²f, P=2N sin²f f'^2, T=sin^4 f/r², F=N r² phi'^2/2 and V=EV+interaction+polarization, energy density e=A+2B+P+T+F+V. Radial pressure pr=2A+2P+2F-e, tangential pt=2B+P+2T-e. Komar mass integrates e^delta(e+pr+2pt)/(3pi). Scalar gradient alone makes zero direct contribution to this combination, but contributes to total mass and changes equilibrium and stresses; not a statement that scalar field energy does not gravitate. Exterior scalar-only Komar volume integrand zero; Komar at boundary differs from ADM by residual f wall stress. ADM mass from asymptotic m, independent asymptotic lapse readout from -g_tt=1-2GM/r+..., and separately integrated mass ledger. Secondary source Q=integral e^delta r² rho(q-alpha phi)/(3pi); verify against asymptotic flux. Flat translational kinetic mass from common action, not a dynamical kick measurement.

## Results

| C=.178396, G=.03, R24 | secondary Q | ADM mass | Komar mass |
|---|---:|---:|---:|
| q=0 |0|1.233413|1.233413|
| q=.25, alpha=0 |.0162936|1.232057|1.232057|
| q=-.25, alpha=0 |-.0162936|1.232057|1.232057|
| q=.25, alpha=1 |.00981192|1.232624|1.232624|
| q=.25, alpha=5 |.00410047|1.233114|1.233114|
| q=.25, alpha=-.25 |.0197637|1.231755|1.231755|

At strong screening Q reduces ~74.8% relative to alpha=0, ADM mass increases ~.086% (less secondary binding). Opposite q gives identical ADM and reversed Q; q=0 has positive mass. C=0 cases show same pattern. For massive q=.25 alpha=0, G=.03, E2=.552780, E4=.656971, EV=.0237547, field=.00138421 including exterior tail .0000522076, interaction=-.00283309, polarization zero; integrated mass 1.23205655444 vs ADM 1.23205655489. Core-only 1.233505, misses ~.118%; binding energy cannot be omitted. G=0 ADM (flat total energy) 1.26088997214, G=.03 ADM 1.23205655489 (2.29% gravitational binding and profile change). Neutral G=0 1.262506→1.233413, similar. These are dimensionless model values, not Newton G or a physical ZPF strength.

## Checks and interpretation

- Max relative mass-ledger discrepancy 4.48e-10 across converged cases; lapse extrapolation discrepancy 1.60e-12 (shares solution, not independent replication). Secondary flux absolute error 7.20e-11, phi infinity residual ≤3.04e-18.
- Komar relative mismatch ≤2.23e-7 for C>0. C=0 R24 discrepancy ≤2.82e-4; R48 ≤3.13e-5. For C=0 q=.25 alpha=0 G=.03 absolute Komar−ADM 2.68543e-4 at R24, 3.35044e-5 at R48, factor ~8; boundary pressure prediction matches to ~1e-10. Finite wall stress, not composition dependence. Two failed tighter refinements do not invalidate standard and wide controls but are unresolved as tighter tolerance tests.
- Flat translation kinetic mass agrees with flat total energy to 9.37e-5 relative (max); no boosted dynamics simulated. For finite G, asymptotically flat conserved energy/momentum relation is conditional on GR, not independently tested with force or acceleration. No full proof of equivalence principle.

**Conclusion:** separating the signed/polarizable secondary source from the universal metric source solves the *accounting* problem under an assumed common action. All contributions, including interaction and self-binding, belong in mass. It does not solve the microscopic origin of metric coupling. GC28D's induced metric kinetic term is compatible with this picture but does not establish why both sectors share a metric. Next question: why general covariance (or an effective universal coupling to conserved stress) arises from the underlying vacuum, rather than inserting it. PI requested that a natural particle-structure cutoff be explored later; no such test was started here. A finite particle form factor may regulate its local fluctuation response without regulating the bulk induced-gravity coefficient, so distinguish those two questions when resumed.
