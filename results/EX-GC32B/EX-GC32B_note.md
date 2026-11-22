# EX-GC32B — Operational geometry from quantum vacuum response

2026-11-22 | ChatGPT | EXPLORATORY | Complete; review/PI closure pending. **PI CLOSED 2026-11-22** (Richard: "I sign off on GC31A–GC33A."). Claude review: review/EX-GC31A-33A_claude_review.md.
[Plan](../../spec/EX-GC32B_plan.md), [source](../../src/exgc32b_response.py), [metrics](metrics.json). Registered plan and LOG before initial production; same-density kinetic-energy control appended to plan after initial gates passed, before those cases. CPU initial 0.518 s, final combined 0.764 s, <2 CPU seconds total, single threaded (NumPy 2.3.5/SciPy 1.17.0). No parent changes; no subsequent test run. Reproduce: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python src/exgc32b_response.py. Archive uses column names and row arrays for compactness; executable emits named rows. The followup repeats 10 initial plane1 records for comparisons, so 42 harmonic + 90 distinct quantum cases, three coherent controls; 40 followup rows including 30 additional cases.

## 1. Question and disposition

Can one distinguish operational geometry, multiple propagation laws and energy-dependent response without fitting geometry to each species? **Methodology gates passed**, not a derivation of geometry or gravity. Stage 1 is a blind-to-adjacency benchmark on known constructed graphs; Stage 2 is an exact finite, nonrelativistic two-sector quantum interaction benchmark with input locality. Neither is a native particle or metric model. Especially significant: varying kinetic energy threefold at fixed source density and fixed H does not change the first response moment, but does change later-time response. This delineates a potential route toward energy-sensitive vacuum dynamics while ruling out equating a density-induced frequency shift with a gravitational metric.

## 2. Harmonic calibration: 42/42

H=1/2(p^T A p+q^T B q), A=M tensor I, B=M tensor L, M=[[1,b],[b,1]], L=.01 I+graph Laplacian. Independent mass-normalized spectral kernel R(t)=F sin(omega t)/omega F^T where F=A^(1/2)V and V diagonalizes A^(1/2)BA^(1/2). Labels randomized (two seeds); pairing of canonical species at a site, clock and units supplied, edges/positions not supplied to estimate_AB. Given R(h),R(2h), A_h=(8R(h)-R(2h))/(6h), C_h=(R(2h)-2R(h))/h^3 -> R'''(0). Richardson at .01,.005; B=-A^-1 C A^-1. Recover edge support by threshold on offsite 2x2 blocks of B, freeze graph before propagation spectrum analysis. Exchange parity known by design. Compare graph-Laplacian eigenvalues with gap-subtracted squared frequencies from inferred D=AB, one independent clock, no per-sector metric.

| Gate | Outcome |
|---|---|
| Edge and all-pairs shortest path recovery | 42/42 exact, including both blinded relabelings |
| A relative Frobenius reconstruction | max 1.05e-8 |
| B relative Frobenius reconstruction | max 4.24e-8 (coarse .02/.01: 6.79e-7) |
| Approximate dispersion speed ratio vs exact | max relative error 1.25e-8 |
| Common cone at 1% tolerance | 14/14 b=0 pass; 28/28 b=.25/.5 correctly fail |
| Static covariance conditional graph recovery | 42/42, BUT cannot distinguish the cone ratio |
| Same-state / energy control | Cq one-component reconstruction vs registered analytic Cq max 2.51e-13 relative; E0 absolute spread max 2.84e-14 within each graph family; full Gaussian state equality analytically in GC32A |
| Coherent excess energy 0,1.005,4.02 | differential impulse response error max 1.43e-16 (fixed quadratic model) |

b=0 ratios 1, b=.25 ratios 1.666666658–1.666666664, b=.5 2.999999963–2.999999987. The full harmonic lattice dispersion has the registered analytical form; the fit uses all nonzero Laplacian eigenvalues (exact for the supplied nearest-neighbor family), not merely first arrival. The continuum gapless coefficient of sqrt(lambda) describes long-wavelength characteristics only; finite, massive, discrete systems have neither a sharp relativistic cone nor an independently established physical metric.

Static Cq is dense (every offsite pair exceeds 1e-4 relative threshold in these graphs). Using the **known harmonic structural prior**, (2Cq)^-2 gives L and the same graph. This is an important qualification: static correlation data can recover adjacency in this special model. Dynamic information is needed to discriminate identical states with different propagation laws, not to prove covariance useless.

Graph-scale diagnostics (mean ball size includes center):
| Graph | r=1 | r=2 | Laplacian gap |
|---|---:|---:|---:|
| cycle12 | 3 | 5 | .267949 |
| cycle24 | 3 | 5 | .068148 |
| cycle48 | 3 | 5 | .017110 |
| torus4 | 5 | 11 | 2 |
| torus6 | 5 | 13 | 1 |
| torus8 | 5 | 13 | .585786 |
| irregular24 | 3.25 | 6.25 | .173893 |

Fit lambda2 proportional to N^(-2/d): d_eff=1.0078 cycles and 2.2631 square tori (three small sizes, finite-size bias, no precision criterion or claim of physical dimension emergence). Graph count and length are input calibration units, not independently derived spacetime scale.

## 3. Nonlinear quantum response: separate onsite shifts, propagation weights and memory

Two otherwise identical probe flavors with bare one-probe H_p=I+L on a supplied 12-site ring. Source bosons H_s=N_s+.5 b^dag L b for mobile cases, or J_s=0 static cases. In the vacuum-to-one-probe response the quartic interaction is exactly representable in a finite 12 or 144 dimensional block. Source states with zero probes are independently energy eigenstates; no metric or gravitational force in H.

- Onsite arm: H_int=g sum_i n_s(i)n_p(i), gives a local density-dependent probe energy (not a link).
- Bond arm (positive method control): H_int=g sum_i n_s(i) W_i, W_i=sum over incident edges (a_i^dag-a_j^dag)(a_i-a_j). This explicitly puts a density dependence in hopping; no geometric emergence claim.

Let Z be the embedding |psi_s> tensor |probe j>. With initial source eigenenergy E_s, Q(t)=Z^dag exp[-i(H-E_s)t] Z is the retarded addition amplitude (probe initially empty). Spectral theorem ensures unitarity of full evolution; projected Q can be nonunitary due to entanglement with source. T=iQ'(0)=Z^dag(H-E_s)Z, S=-Q''(0)=Z^dag(H-E_s)^2 Z. Independently recover T,S via finite time samples at .01,.005 and Richardson. For both arms T=H_p+g sum_i rho_s(i)W_i with W_i=|i><i| onsite. V=S-T^2=Z^dag(H-E_s)(I-ZZ^dag)(H-E_s)Z >=0 is the second-moment excess, a diagnostic of coupling beyond one time-independent projected generator, not a metric. This exact identity follows from Z^dag Z=I (orthonormal probe basis). T's offdiagonal elements infer edge weights; separate residual diagonal after subtracting common gap and weighted graph Laplacian. Full Q(t) and V provide independent information; time evolution by T alone is not generally exact.

Stage 2 initial 60 rows all pass: max T error 1.66e-7, S error 7.41e-8, direct independent projected-moment identity relative error 2.84e-16, mobile source energy eigenresidual 4.69e-16. All g=0 time responses to bare baseline agree within 7.0e-15. Stable nonnegative couplings, Hermitian time-independent full H; total energy conserved during evolution; E_s computed from no-probe source H (not inferred from geometry). Static localized n=2 is exact conditional invariant block; neither fixed source nor prescribed ring is a dynamical emergent particle.

Illustrative results at g=1 (unit baseline edge weight; Q difference = relative Frobenius norm vs bare t=2; values rounded):
| Source | E_s | Onsite: peak onsite shift | Onsite Q difference | Bond: peak edge shift | Bond Q difference |
|---|---:|---:|---:|---:|---:|
| Vacuum | 0 | 0 | 0 | 0 | 0 |
| Local 1 | 1 | 1 | .3832 | 1 | .5696 |
| Local 2 | 2 | 2 | .6527 | 2 | .7777 |
| Uniform 1 static superposition | 1 | .08333 | .1624 | .16667 | .3347 |
| Plane wave m=1 | 1.13397 | .08333 | .1532 | .16667 | .3386 |
| Standing wave m=±1 | 1.13397 | .16667 | .1730 | .29167 | .3389 |

Onsite has no inferred edge change (max spurious 3.96e-8); bond has no residual onsite shift (max 2.09e-11); max spurious nonedge weight 4.55e-7 due to finite-difference truncation, far local-source edge shift <1.52e-8. Bond local1 edge increment .25/.125 or 1/.5 gives 2:1 (within <3e-7 absolute derivative error); same bare probes can see different effective link weights simply by interaction strength. Equal g produces equal response by construction, not derived universality. The local changes vanish outside incident edges at first moment: not an inverse-distance gravitational field. Q norm is not a monotonically calibrated force; e.g. local bond at t=2 is .6093 for g=.5 vs .5696 for g=1. Equal energy localized vs uniform or plane vs standing may have different spatial distribution/stress; differing local response is expected and not alone an equivalence-principle violation.

## 4. Registered followup: total energy vs density on the SAME Hamiltonian

This was added prospectively after the initial methodology gates passed. Fix J_s=.5, n=1, choose ring Bloch eigenstates m=0,1,3,6. All have identical site probabilities 1/12 in exactly the same microscopic H. E_s=1+.5(2-2 cos(2pi m/12)), so E_s=1,1.1339746,2,3. This changes kinetic energy without changing density or H; momentum/stress changes as well, so not a test of full relativistic equivalence.

In all 40 comparisons (30 new plus 10 original plane rows), first moment relative difference from k=0 is <=1.53e-9 (<1e-5 criterion); T and S error max 2.44e-8 and 1.13e-8. Exactly, T=H_p+g sum W_i/12, independent of kinetic energy. V is also independent here for the plane states: since (H_s-E_s)|psi>=0, V is covariance of the interaction by itself, so density probabilities suffice. Higher time response does distinguish the states, however. At g=1, Q(t=2) relative difference from m=0:
| m | E_s | Onsite full-response difference | Bond full-response difference |
|---:|---:|---:|---:|
| 0 | 1 | 0 | 0 |
| 1 | 1.13397 | .01587 | .11575 |
| 3 | 2 | .04275 | .25844 |
| 6 | 3 | .05246 | .29745 |

The exact V Frobenius norm is .264619 onsite and 3.00617 bond for all four momenta at g=1. Thus even *two* short-time moments do not reveal all kinetic-energy-dependent dynamics. Full temporal kernel changes (source mobility, phases and interaction dynamics). Cannot say full response is sourced simply by energy, because momentum/internal state also change; no shared cone in this nonrelativistic benchmark, no active gravitational mass or force law. In particular the first-moment density-dependent bond rule misses kinetic energy, so its apparent common-link deformation is not a general energy-sourcing mechanism. This is arguably the most informative negative finding.

## 5. Interpretation and next direction

The GC32A hypothesis is sharpened:
- A quantum state is not a sufficient dynamical description of a common operational metric.
- Retarded response can reconstruct a local generator in controlled circumstances, but chosen observable basis, clock, and Hamiltonian restrictions matter.
- A change in response need not imply spatial geometry. Decompose onsite vs link changes and frequency dependence; even an explicitly constructed common link shift misses total kinetic energy here.
- Total energy has to be defined independently. The source-motion sensitivity of the full kernel is encouraging for a dynamical mechanism, but does not establish that the kernel obeys an energy/stress-based law, that the response is universal, or that gravity follows.

**Recommendation (proposal only): GC32C — vacuum-selected relational geometry.** Rather than enlarge this supplied ring benchmark, evaluate a small permutation-invariant dynamical-link substrate with no fixed spatial adjacency, linking to the Sakharov line. For example, allow undirected bond variables on all abstract label pairs, with local edge degrees as dynamical variables; couple vacuum oscillator modes through the resulting relational Laplacian. Compute the mode ground energy independently and ask whether vacuum fluctuation energy differences favor spatially extended connected low-dimensional networks versus random, clustered or disconnected networks. Include bare link-only or constrained edge-count controls (constraints are inputs, not emergent laws), avoid prescribing a target lattice. Pure positive oscillator zero-point energy is monotone with edge addition for positive stiffness; without additional stabilization it favors empty links, so test rearrangement at fixed edge budget first, not assume a stable vacuum. A tentative next discriminant is whether changing independently tracked excitation energy and internal state changes link equilibrium/retarded response through the same energy functional, not a hand-inserted source term. Need to distinguish kinetically trapped patterns from stable low-dimensional vacuum and establish perturbation robustness before particle tests. No reason yet to claim this construction is novel or will succeed; this is a concrete bounded screen. If no vacuum-driven preference for useful geometry, stop rather than tune a lattice into existence.

**Disposition:** complete, all numerical methodology gates met; no gravity success claim. No GC32C or 31B2/B3 authorized or run. Review and PI closure pending.
