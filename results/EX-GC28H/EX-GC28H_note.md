# EX-GC28H — even vacuum coupling and fifth-force suppression (EXPLORATORY)

**Summary:** A common fractional mass-squared coupling to a loop-induced secondary scalar removes the odd-source cancellation in 28C but yields different regulated source/energy ratios for unequal mass species. The collective scalar still needs mass tuning, and unlike the symmetric ±s model it also needs tadpole and cubic tuning for a stable massless stationary point. Strong intrinsic polarization can suppress the secondary monopole by orders of magnitude in a fixed-core Skyrmion, but finite polarization does not enforce exactly zero or proportional source. A small secondary mass screens the force at sufficient separation. Neither mechanism establishes universal gravity; the arbitrary G used in 28G cannot be extrapolated to actual gravitational strength.

Design/run ChatGPT; independent review and PI closure pending. PI authorized 28H only; natural cutoff 29A deferred. Plan `spec/EX-GC28H_plan.md` and LOG committed before run. Source `python src/exgc28h.py`; 24 loop/detuning rows, 32 source rows, 10 fifth-force baseline rows and 80 bound rows, 150 fixed-core screening rows. 1.68 CPU seconds, single thread, no deviations. Reuses 28C proper-time functions, 28E neutral profile solver and committed 28G records. All new numerical rows retained.

## 1. Even loop coupling

Microscopic mass squares a_i(1+s), a=(1,1) or (1,4), no bare kinetic term. T=sum a_i I(a_i)/2 (one-loop tadpole), ucrit=sum a_i² B_i(0)/2; Z=sum a_i² B'_diff(0)/2. A linear local counterterm must cancel T, a quadratic one must tune K(0)=0; without quadratic counterterm the reference curvature is negative. Z positive. Cubic V'''=sum a_i³ exp(-a_i/Lambda²)/(32pi² a_i) is positive, verified by central second derivative of tadpole at ±.001 to 3.95e-7 relative. At tuned zero curvature and nonzero cubic the stationary point is not a local minimum; a separate cubic cancellation (and positive quartic) would be required. No stable massless branch claimed. Detuning by Z mu² gives positive low-momentum quadratic pole mass close to mu (records). Common mass coupling is dilaton-like but not the full universal conformal metric coupling.

| a=(1,4), Lambda | tadpole | ucrit | Z | cubic |
|---:|---:|---:|---:|---:|
| 2 |.0140800|.0144206|.00118753|.0211029|
| 5 |.265987|.0798799|.00230578|.0462122|
| 10 |1.38171|.148620|.00255055|.0518089|
| 20 |6.05729|.221715|.00261625|.0533149|

Defects with infinitesimal nonnegative lambda_i shift mass squares; E=sum lambda_i I_i/2, J=sum a_i lambda_i B_i(0)/2. At Lambda10 the pure light vs heavy source/energy .0425193 vs .125654 (ratio 2.955); both positive. Equal a=(1,1) gives equal ratios by construction. At Lambda2 values .504260/1.47738; Lambda20 .0137626/.0425193. Loop source is regulator dependent; Born energy is quadratically divergent, so these are not physical total mass or equivalence-principle predictions. Even sign symmetry removes the particular odd cancellation, not composition dependence. This narrows, not extends, 28C's original conclusion; agree with Claude's critique.

## 2. Fifth-force accounting correction to 28G

In radial normalization E=(1/3pi) integral density dr and Z=1, scalar tail phi=3pi Q/r, pair energy -3pi Q1 Q2/(Z r), vs metric -G M1 M2/r (weak long-distance). Like sign attracts, opposite repels. Identical-pair relative force 3pi (Q/M)²/(ZG); at 28G G=.03 ratios unscreened .0549444, strong-screen .00347, anti-screen .0809. Neutral has zero. This is an extra non-universal force, even when ADM accounting is correct. The relative strength scales 1/G; no physical conversion from dimensionless G. For illustrative required ratio epsilon, |Q/M| <= sqrt(epsilon ZG/(3pi)): at G=.03 epsilon=1e-13 bound 1.78412e-8; at G=1e-12 bound 1.03006e-13. Neither is a claimed experimental limit. For finite range mu^-1, pair force multiplied by (1+mu r)exp(-mu r) in point/tail approximation; using 28G unscreened scalar amplitude fixed, require mu r=30.48 to reach 1e-13 at G=.03, 55.19 at G=1e-12. Actual massive scalar amplitude and dressed energy must be recalculated, not held fixed for physical prediction.

## 3. Fixed-core high-polarization test

Neutral radial Skyrmion profile f held fixed (C=0 or .178396), rho=(1-cos f)^2. Solve u=r phi, -u''+(mu²+alpha rho)u=q r rho. Positive operator alpha>=0, u(0)=0, exterior Robin u'(R)+mu u(R)=0. Symmetric tridiagonal variational discretization with endpoint half weight; finite h and f(R)=0. Exterior energy included, QY=u(R)exp(mu R)/(3pi). For mu>0, QY is Yukawa amplitude, not a conserved Gauss source. Compare independent integral moment using sinh(mu r)/mu and exact discrete summation-by-parts. Max relative continuum moment discrepancy 3.82e-5, absolute discrete identity error ≤4.61e-13, field/interaction ledger ≤8.44e-15. No high-alpha nonlinear profile relaxation (fixed core only). q=.5 doubles field/source at same alpha; energy shifts quadratically.

| alpha, C=.178396 q=.25 mu=0 R24 h=.05 | Q/E | force ratio at G=.03 | total energy |
|---:|---:|---:|---:|
| 0 |.0143353|.0645600|1.260912|
| 1 |.00853114|.0228646|1.261592|
| 5 |.00351552|.00388266|1.262165|
| 25 |.00104621|.000343863|1.262422|
| 100 |.000334449|.0000351406|1.262484|
| 1000 |.0000467355|6.86188e-7|1.262504|
| 10000 |.00000624993|1.22716e-8|1.262506|
| 1000000 |.000000102261|3.28525e-12|1.262506|

Very strong polarization suppresses Q/E ~1.40e5 times but still misses illustrative 1e-13 force bound even at large toy G=.03. At G=1e-12 the largest-alpha massive-core q=.25 ratio is ~.0986 (not tiny). At q=.5/C=0 alpha=1e6 ratio at G=.03 =3.81e-11, illustrating failure of universal suppression. Positive finite alpha and nonzero q give positive u and tail, not exact cancellation (positive Green function); alpha→infinity can drive amplitude toward zero, but also drives scalar dressing energy toward zero at fixed q. At alpha=1e6 C>0 q=.25 Q=1.29105e-7 vs 1.29106e-7 at h=.025, R48 1.29105e-7. C=0 longer tail: R24 Q=2.14437e-7 vs R48 2.31131e-7 (~7.8%); high alpha samples more of the outer profile, so domain effect explicitly unresolved. No extrapolation beyond tested alpha. For C>0 alpha=1e6, q=.25, mu=.05/.2 QY=1.64556e-7/3.38429e-7, larger Yukawa amplitude than massless but extra exponential suppression at distances large compared to range and particle size. Not equivalent to exact zero source or universal screening at all distances.

## Decision

Accept Claude's 28B–G review: correct 28B wording energy 1/d³, force 1/d⁴; 28C odd source built into coupling and Goldstone angle derivative sourcing is expected under its symmetry; even coupling tested here and also fails to establish universal physical sourcing. 28G metric mass accounting remains sound but its massless secondary scalar is a distinct fifth force. A finite-range field or a symmetry forbidding its monopole could keep local dressing while avoiding long-range composition-dependent forces; neither arises automatically from screening. No further numerical precision chase. Natural cutoff EX-GC29A is next recommendation (not run here); distinguish particle form factors/spectral responses from bulk vacuum induced EH stiffness. Common metric emergence remains open.
