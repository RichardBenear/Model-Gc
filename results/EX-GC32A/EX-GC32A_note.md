# EX-GC32A — Vacuum dynamics, operational geometry and gravity

2026-11-22 | ChatGPT | EXPLORATORY | paper-only | complete; independent review and PI closure pending.

## 1. Question and scope

Can physical geometry be an output of vacuum dynamics, so that particle energy and gravitational response are properties of the same substrate? This is the next branch after the GC30 common-metric screens and GC31 native-rewrite studies. See [registered plan](../../spec/EX-GC32A_plan.md). No production simulation, negligible analytical compute; within 2 core-h. No changes to parent ModelG, no GC31B2/B3 or GC32B run. The PI's interest in Haramein and Heim is included without assuming their claims.

The prior results distinguish three tasks: (i) derive geometry, (ii) explain universal operational propagation on it, (iii) explain why conserved energy sources changes of it. GC28D already obtained a conditional induced Einstein–Hilbert coefficient on a supplied metric; GC28G confirmed the distinction between signed secondary scalar charge and positive metric energy in a model with metric coupling supplied. GC29A separated spatial form-factor rolloff from a universal vacuum UV cutoff. GC30 showed topology, microscopic exchange invariance and weak finite-scale RG attraction alone do not ensure a common metric. GC31 showed that schedule invariance, shared structural reach and active relational turnover alone do not do so either.

## 2. Prospectively registered criteria and literature screen

G1: derive geometry (order, scale, dimension and clock) rather than give it to the matter model.
G2: native persistent localized excitation, not a prescribed defect.
G3: independent conserved energy, not a quantity defined by gravitational response.
G4: derived feedback of energy onto geometry.
G5: common effective principal cone and response across physically distinct sectors.
G6: finite viable vacuum stiffness/UV mechanism and massless tensor coupling (not simply a long-range scalar).

P = demonstrated within the stated model; C = conditional, assumed, or only part of criterion; U = unestablished by the selected construction; F = failure of a specified implication in the tested model. No aggregate scores: the criteria are not interchangeable, and U does not mean a universal no-go.

| Candidate (selected construction) | G1 | G2 | G3 | G4 | G5 | G6 | Relevant capability / limitation |
|---|---|---|---|---|---|---|---|
| Causal sets (order + counting) | C | U | U | U | C | U | Manifoldlike causal order plus volume can recover continuum geometry up to applicable assumptions; discreteness supplies a counting scale when density specified. Retarded scalar operators and curvature estimators exist. Dynamics selecting a physical causal set, native particle energy, and multi-species probe equivalence are separate. |
| Quantum causal networks / free quantum cellular automata | C | U | C | U | C | U | Discrete unitarity, local update and symmetry can yield a Dirac limit. Neighborhood/Cayley group, periodic steps and isotropy assumptions supply geometry; quasienergy is modulo the step frequency. Distinct excitations need not share the same limit. Not a background-free derivation of matter or gravity. |
| Entanglement geometry (Van Raamsdonk) | C | U | C | C | U | U | Connected bulk geometry linked to entanglement in holographic quantum gravity, under the duality and chosen state. Does not establish a unique physical clock, common metric across arbitrary sectors, or a particle-derived cutoff from an arbitrary state. |
| Entanglement equilibrium (Jacobson) | C | U | C | C | C | C | First-order local vacuum equilibrium at fixed volume yields semiclassical Einstein equation for conformal quantum matter under entropy and geometric assumptions. Uses local Lorentzian causal diamonds and stress-energy. Does not derive these from pregeometry or independently determine G. |
| Sakharov vacuum-induced stiffness | C | U | C | C | C | C | Integrating fluctuations in curved background induces an effective curvature action. Can link vacuum dynamics to metric elasticity; common metric is an input, coefficient and vacuum term regulator/spectrum dependent. No native particles or unique physical cutoff established. |
| Haramein (2013 holographic mass model) | C | U | C | U | U | U | Proposes Planck-scale geometric counting tied to a Schwarzschild mass-radius relation and a proton mass formula. Geometric scaling and numerical agreement are not independent conserved microscopic energy or dynamical stress-energy backreaction. G is used to construct Planck units and the Schwarzschild relation; this is not an independent derivation of G. Proton charge radius is not a demonstrated gravitational horizon. Extra physics would have to identify what gravitates, predict propagation, and reproduce weak-field observations. |
| Heim (1977 particle theory proposal) | C | U | U | U | U | U | Higher-dimensional geometrical carrier and proposed elementary-particle description. The dimension/geometry is postulated; a proposed spectrum is not an operational derivation of energy conservation, universal cones and tensor response. No endorsement or general disproof; a reproducible dynamical model is needed for these tests. |

These are judgments about what the cited construction establishes for *our* gates, not rankings of the entire research programs. In particular, entanglement geometry may specify additional dynamics through a duality, and causal sets might acquire geometry-selecting quantum dynamics in other models.

Weinberg–Witten constraint: a Lorentz-invariant theory with a Lorentz-covariant conserved stress tensor cannot contain massless spin > 1 states carrying its nonzero energy-momentum. An emergent-gravity proposal must identify which premise does not hold, and still recover observed low-energy Lorentz invariance. A fundamental discrete order or non-Lorentz-invariant quantum network could avoid the microscopic Lorentz premise; a genuinely gravitational/holographic description need not have a local gauge-invariant covariant gravitational stress tensor. Neither escape establishes a graviton. Sakharov/Jacobson start with geometric degrees of freedom or semiclassical GR rather than construct a composite graviton in an ordinary flat-space QFT. Haramein/Heim have no demonstrated relevant spin-2 emergent mechanism on this screen; assigning an escape without a precise dynamical derivation would be premature.

## 3. Exact analytical arm — same vacuum, different operational cones

Use a finite periodic demonstration line with N sites, real canonical coordinates q and p, two components per site, hbar=1, site spacing and clock units fixed. L=m^2 I+graph Laplacian, m>0 for a unique normalizable ground state. Let

    H_b = 1/2 (p^T A p + q^T B q),
    M_b = [[1,b],[b,1]],  |b|<1,
    A=M_b tensor I_N, B=M_b tensor L,
    K=I_2 tensor sqrt(L).

Both A and B are positive; B=KAK. This has genuine off-diagonal microscopic bilinear kinetic and potential couplings at b≠0, exact exchange symmetry, and one shared time. It is Gaussian, hence diagonalizable and not a model of nonlinear scattering. The graph, dimension, spacing and time have been deliberately supplied; the goal is a counterexample to the *sufficiency* of equal-time vacuum data, not a full pregeometric model.

For psi(q) proportional to exp(-q^T K q/2), p^T A p psi = [tr(AK)-q^T KAK q] psi. The quadratic term cancels the potential, giving H_b psi = tr(AK) psi/2. Stability implies this Gaussian is the unique ground state. Its full covariance is

    <qq^T> = K^(-1)/2,
    <pp^T> = K/2,
    <{q,p^T}>/2 = 0.

This is the **same entire state for all b**, not just one density statistic. All equal-time correlations, and all entanglement measures for any fixed subsystem partition, are identical. Even the unshifted total vacuum energy is E0=tr(sqrt(L)), independent of b since tr(M_b)=2. This includes the correlations between all sites. A specification of this Hamiltonian family or additional dynamics can distinguish the models; the ground state alone cannot.

The symmetric/antisymmetric normal modes have M eigenvalues 1±b. Fourier-transforming the given lattice yields

    omega_±(k) = (1±b) sqrt(m^2 + 4 sin^2(k/2)).

Take the longwave gapless limit m→0 (remove the zero mode, or keep m<<|k|<<1 on a large lattice). The low-energy characteristics are c_±=1±b in the common clock units. Do not confuse these with massive group velocity or claim a sharp microscopic front in the continuous-time harmonic lattice.

| b | c+ | c- | ratio | ground state | total vacuum energy |
|---:|---:|---:|---:|---|---|
| 0 | 1 | 1 | 1 | same | same |
| 1/4 | 1.25 | 0.75 | 5/3 | same | same |
| 1/2 | 1.5 | 0.5 | 3 | same | same |

The ratio is invariant to common changes in clock and length units. No single low-energy principal light cone exists for the two modes when b≠0. One cannot select a separate time per species and call the result a common geometry. This sharpens GC30A's exchange-symmetry example and GC31's common structural-front counterexample: even exact quantum state and total vacuum energy agreement do not enforce speed agreement. It does not exclude holographic or other entanglement-based constructions that supply extra dynamical constraints. Model b=0 is a deliberately imposed common-cone positive control.

### 3.1 Retarded response as independent data

A delta-function force coupled to q at t=0 shifts p, so the impulse response R(t) is the matrix mapping the initial momentum shift to q(t). Equivalently R(t)=i theta(t) [q(t),q(0)^T]. With qdot=A p and qddot=-AB q,

    R(t) = (AB)^(-1/2) sin[(AB)^(1/2)t] A,  t>0,
    R(t) = t A - t^3 A B A/6 + O(t^5),
    R'(0+) = A,
    R'''(0+) = -A B A,
    D = - R'''(0+) [R'(0+)]^(-1) = AB,
    B = - A^(-1) R'''(0+) A^(-1).

The Gaussian family has the same state but different retarded kernels. Thus response, with fixed observable basis and independent common clock, supplies information absent in the state. It determines a generator here, not yet a physical metric. Its spatial/spectral structure can be tested using blind reconstruction in a known benchmark. Labels, canonical normalization and the clock supplied to the example are assumptions; a more ambitious candidate must account for them.

### 3.2 Energy without geometric feedback — exact negative control

For any state of this fixed quadratic Hamiltonian, q(t) is linear in q(0),p(0). Their canonical commutators are c-numbers, so the retarded response R is state independent. A coherent excitation can deposit arbitrarily more positive conserved energy (relative to E0), and exchange symmetry can be kept by exciting the symmetric mode; R and its principal cones still do not change. No energetic backreaction on operational geometry occurs. This exact consequence also holds for other excited states in the same Hamiltonian, not just coherent displacements.

This is not a statement about a dynamical geometry with a quadratic *matter* sector and independently coupled gravity. It is a failure of the fixed Gaussian substrate itself to produce energy-dependent response. A candidate needs something beyond this: dynamical relational variables or genuinely nonlinear/interacting response, or another independently derived mechanism. Nonlinearity is not by itself universality, nor is changing a susceptibility necessarily curvature.

### 3.3 Registered criteria: analytic model

- G1 C: given graph time/dimension/length; response identifies generator conditional on those inputs; not emitted physical geometry.
- G2 U: free modes and coherent packets, no demonstrated persistent native topological particle.
- G3 P (within finite model): time-independent stable H conserved; excitation energies omega_± computed independently.
- G4 F for fixed quadratic substrate: increasing energy cannot change commutator response.
- G5 F for b≠0: exact equal-state/exchange-symmetric two-cone witness; C for b=0 imposed control.
- G6 U: graph regulator given, no tensor dynamics or derivation of physical stiffness.

No test candidate passes the full gate. The result is constructive: distinguish correlations (state), response (dynamics), and energy-induced response change (candidate backreaction). Paper-only analytical arm executed as registered.

## 4. Consequence for the Sakharov hypothesis

A useful possible combination, not a derivation, is: relational quantum vacuum dynamics -> common operational metric -> induced stiffness under changes of that metric -> conserved energy as source. This is consistent as a research direction with Jacobson's conditional entanglement-equilibrium argument and Sakharov's induced-action proposal. It requires deriving the arrows rather than combining assumptions from the literature.

In four dimensions a schematic regulated effective action on an *already supplied* metric contains terms of scale Lambda^4 for the volume and Lambda^2 R for curvature, plus higher derivatives/logarithms. The numerical coefficient and sometimes its sign depend on field content, couplings, regulator and conventions (see GC28D); a physical cutoff and small renormalized cosmological term are independent requirements. A frequency cutoff tied to a preferred inertial frame is not automatically Lorentz compatible. The GC29A particle form factor is not that vacuum regulator. One must demonstrate viable tensor degrees of freedom, positive stiffness, universal coupling to total conserved energy including binding, and a weak-field limit with the correct long-range relation. The static vacuum energy itself should not simply be identified with the gravitational source. The Gaussian example shows precisely why static energy/density cannot determine the dynamical response even when equal.

## 5. Recommendation — GC32B proposal only (not authorized/run)

First perform cheap operational response tomography in small quantum oscillator networks with blinded site adjacency (registered controls), or a finite exactly tractable quantum network. The observer receives canonical observable labels without coordinates/neighbors; all probes use one independent common clock. Register a deterministic inference rule before observing outcomes:

1. Recover the response generator from retarded impulse data, compare independently with the state covariance graph/spectral estimates. Blind graph labels and require relabel covariance; infer candidate distances and dimension from connectivity/spectral scaling, not by redefining lengths to match each probe's velocity. Graph recovery is a methodological control, not emergence if an input Hamiltonian already has locality.
2. Positive imposed common-kernel control; same-state/different-H negative control above; independently prepared unequal-cone probes, and matched-energy excited-state control in quadratic H. Require both mode dispersion and principal characteristics measured against one geometry; distinguish reconstruction from selection of a common geometry. Finite network limitations explicit (no sharp spacetime cone or continuum at small N).
3. If these controls work, select a small nonlinear or dynamical relational candidate where the retarded operator may depend on the excitation, preferably without an assumed spatial neighbor graph. Compare configurations with equal conserved excess energy but distinct internal structure and controls at unequal energy. Measure change in propagation / inferred generator independently of the excitation energy (do not fit geometry to force). Demand reproducible common response across probe sectors, check localization, energy ledger, stability, far-field decay and sign where appropriate. Reject species-dependent index changes as evidence of GR. The independent clock provision is a calibration device, not an emergent time claim. This is a second gate, not automatic continuation to an expensive run.

An all-to-all permutation-invariant relational Hamiltonian with dynamical links may be an interesting candidate after the benchmark; spontaneous selection of locality/dimension and robust universal propagation are not guaranteed. No fixed rewrite gadget or prescribed geometric metric should be added just to obtain a positive score. Only after state-dependent geometry can be shown in a common response should we ask whether its effective stiffness is induced by the vacuum spectrum, whether the source scales with conserved energy, and whether the response has the tensor and long-range properties of gravity.

## References (primary works)

- Bombelli, Lee, Meyer, Sorkin (1987), [Space-time as a causal set](https://doi.org/10.1103/PhysRevLett.59.521).
- Benincasa & Dowker (2010), [The scalar curvature of a causal set](https://arxiv.org/abs/1001.2725).
- Johnston (2010), [Quantum fields on causal sets](https://arxiv.org/abs/1010.5514).
- D'Ariano & Perinotti (2014), [Derivation of the Dirac equation from principles of information processing](https://doi.org/10.1103/PhysRevA.90.062106).
- Van Raamsdonk (2010), [Building up spacetime with quantum entanglement](https://arxiv.org/abs/1005.3035).
- Jacobson (1995), [Thermodynamics of spacetime: the Einstein equation of state](https://doi.org/10.1103/PhysRevLett.75.1260); (2016), [Entanglement equilibrium and the Einstein equation](https://doi.org/10.1103/PhysRevLett.116.201101).
- Sakharov (1967), [Vacuum quantum fluctuations in curved space and the theory of gravitation](https://www.mathnet.ru/eng/dan33444).
- Haramein (2013), [Quantum gravity and the holographic mass](https://journalpsij.com/index.php/PSIJ/article/view/158), Physical Review & Research International 3(4), 270–292 (2013).
- Heim (1977), [Vorschlag eines Weges zur einheitlichen Beschreibung der Elementarteilchen](https://doi.org/10.1515/zna-1977-3-404).
- Weinberg & Witten (1980), [Limits on massless particles](https://doi.org/10.1016/0370-2693(80)90212-9).

**Disposition:** paper test complete, review/PI decision pending. GC32B is a proposal, not an executed test. No numerical precision campaign; no claim to have derived physical geometry, inertia, G, or gravity.
