# EX-GC24A — reactive medium response versus head-on drag (exploratory)

PI approved revised tests, following Claude review `5fd7aa5`. Base proposal `c2b3ce6`. Protocol `spec/EX-GC24A_plan.md`; run `python src/exgc24_impedance.py` from repo root (A5 diagnostic needs EX-GC23 saved rest arrays). Machine-readable output `diagnostics.json` (14 spectral configurations × 6 frequencies, 14 steady controls, 3 independent time integrations, 6 actual 1/G harmonic evaluations, 4 centroid diagnostics). Two complete runs total ~25 CPU s, plus small reporting overhead. A1–A4 complete; A5 instantaneous diagnostic complete, dynamic refinement deferred; B stopped by preregistered conditional gate. No free particle with derived inertia was evolved.

## Main result

For a prescribed Gaussian sink in the existing telegraph medium, positive low-frequency reactive inertia is reproducible at weak damping; negative reactive response at strong damping is also reproducible. At fixed shape and G=1 in the force proxy, inertia and wake drag scale exactly as S², while G sourcing and the external push scale as S. The existing head-on closure overwhelms the reactive mass: m/(b_bare+b_wake) ≈ .0433 for reference S=5.32 and .00217 at S=.266, far below the 20-unit EX-GC23 tracking interval. **Conditional negative result for the existing rest-frame absorption closure.** The wave equations alone do not derive that closure, and the proxy is not physical momentum conservation.

## A1–A2: analytical/spectral results

For exp(-iωt), Z=-F/V, b=Re Z, m=-Im Z/ω. Let a=c²k², W=2|s_hat|² k_z²/L³. Exact finite-box low-frequency expansion gives b0=Σ W c²/(a² τJ) and m0=Σ W c²/a² [1−1/(a τJ²)]. The damping-squared correction is essential: the ratio is τJ only in the weak-damping limit. The lossless Gaussian continuum answer S²/(6π^(3/2)c²σ)=.2647 at reference settings. Order of limits matters: the correction has enhanced infrared sensitivity.

At ω=.003, S=5.32, σ=.8, c=2, L=32, h=.5:

| τJ | m | b | m/b |
|---:|---:|---:|---:|
| .25 | -4.88762 | .92544 | -5.28139 |
| 2 | .15131 | .11577 | 1.30699 |
| 20 | .23074 | .01158 | 19.93045 |
| 200 | .23153 | .001158 | 199.99058 |

The reviewed Claude τJ=20 and 2 ratios reproduce. For τJ=20, m at ω=.003 is .1986/.2307/.2463 for L=16/32/64. Halving h at L=16 changes the lossless mass by 4.17e-13 absolute. At reference τJ=20, varying S=.266,2.66,5.32,10.64 gives m=.0005768,.05768,.23074,.92294; doubling c quarters lossless m; changing σ changes m approximately inversely with width plus finite-box correction. Neither changing τJ nor refining h supplies the missing scaling. Even the largest conditional m/S in the one-factor sweep (c=1) is .172 time units. No claim over untested parameters or realistic core shapes.

Additional actual recoil linearization (post-protocol extension, no tuning): retaining 1/G in the force and differentiating the moving source, J and denominator gives m=.231242 and b=.0249605 at ω=.003 for S=5.32; at S=.266, m=.000576907 and b=.000030605. Thus the weak-field approximation accurately predicts mass here but can underpredict wake drag by >2× at strong S despite small depletion. The ratio m/b is then ~9.26 rather than 19.93. For the weak sink the ratio is ~18.85. This is distinct from bare drag.

## A3: independent time integration

Evolve G and div J first-order ODEs by RK4 with the full sinusoidally translated Gaussian transform (not the harmonic response formula), starting from static equilibrium. Fit last two of six periods to F=-b v-m a. L=12, h=.75, τJ=2, ω=.5. Spectral b=.1116286380, m=.1678844807. A=.02, dt≈.020 returns b=.1116329782, m=.1678829952. Halving dt hardly changes results. Halving A to .01 gives b=.1116297204, m=.1678841026. Worst fractional discrepancy .00389%; all within 1% preregistered gate. Independence is time integration versus frequency-domain algebra, sharing spatial Fourier discretization, not an independent spatial stencil or physical validation.

## A4: absorption comparison

EX-GC17 measures drag in token model. EX-GC19/20 continuum extrapolates this as explicit -drag*V in the force per consumption, with localization weight s/(s+.001 max s). For prescribed uniform v, difference between enabled/disabled force is identically -S v (or -integral s wgt v). This calculation checks the inherited closure, **not** whether absorbing wave G would physically produce such a term. At v=.05 reference S=5.32, actual G wake force=-.00124897; bare=-.266 (local-weight bare=-.262271). At S=.266 wake=-1.5314e-6, bare=-.0133. +/- velocity and zero controls pass; uniform flux +/- .02 changes external force independently. At positive actual harmonic m, b_total=5.34496, t_coast=.04326; weak t_coast=.002169. Both fail conditional <1 and <1/20 of tracking interval gates. A simple local exponential decay using this proxy is only an illustrative extrapolation: response is dispersive and short times probe outside measured low-frequency band. Do not claim actual decay trajectories. A missing independent wave absorption law remains open; removing the term solely to get coasting is not justified. Stop B under the tested closure, not a universal no-go theorem for all field inertia.

## A5: slip diagnostic and cost

Using original EX-GC23 settled W and V arrays and cubic interpolation, compare instantaneous W centroid velocity increment from uniform boost .05 using native centered transport stencil. Strong K=.2 gives ratios .99999926/.99999903 at h=.25/.125; weak K=.01 .99999867/.99999827. This rules out a simple immediate uniform-boost centroid derivative defect, not accumulated evolution slip. Full independent .125 settling and 20-unit coasts are estimated at 32× the .25 CPU time: ~14 core-hours for strong+weak pair (rest and kick combined). Full dynamic convergence unresolved, deferred as stated in protocol, no claim that the 2.5% slip is fixed. Raw rest arrays previously supplied in EX-GC23_settled_fields.zip.

## Next candidate

EX-GC25 should first specify a reciprocal Lorentz-compatible coupled theory with a stable localized rest configuration and an explicit energy-momentum ledger; test invariance of boosted solutions and E(v), P(v) rather than simply switch off drag. An alternative is internal/core inertia (possibly topological); retain it as a competing hypothesis. Test gravitational source proportionality only after the inertial coefficient is independently derived. A Lorentz-invariant scalar wave equation alone is insufficient when chemistry, damping, consumption and reservoir define a rest frame. No orbital tests yet. No physical mass-energy equivalence or universal free fall established.
