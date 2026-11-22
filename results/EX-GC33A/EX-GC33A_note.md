# EX-GC33A — Quantum vacuum phase, relational geometry and universal energy response

2026-11-22 | ChatGPT | EXPLORATORY | paper screen complete; independent review and PI closure pending.
[Plan](../../spec/EX-GC33A_plan.md). Registered plan/LOG before analytical work. No production simulation; negligible arithmetic, within 2 core-hour cap. No parent changes. GC31B2/B3 remain deferred. No further production authorized.

## 1. Purpose and prior-work disposition

GC32A proved static vacuum state and even total vacuum energy do not fix a common propagation cone. GC32B reconstructed conditional operational graph geometry from response and found later-time response can depend on source motion even when its density and first two response moments agree. GC32C showed that a scalar positive oscillator ground-energy potential favors deletion, clustered networks and disconnected blocks rather than the reference extended lattices; the link force missed uniform massive excitation energy.

The PI notes that Quantum Graphity has already been explored in earlier threads (possibly with Claude). The accessible project-history retrieval did not identify a specific completed Graphity implementation or its review. It did surface a prior StochasticVacuum DN-01 equal-rate connected-graph pivot null (reported featureless with short distances), and CRR's mixed matter/vacuum-rule obstruction. These are pointers, not independently re-audited tests. Respect the PI's preference: do not restart a valence/loop parameter search or presume these experiments are new. Use Graphity as established prior art and ask a different discriminating question.

Central distinction: (i) what dynamical law selects the vacuum geometry; (ii) whether different physical sectors share its characteristic operator; (iii) whether conserved energy/stress affects that geometry; (iv) whether the geometry response is the observed spin-2 gravitational interaction. They are not interchangeable. G1 geometry without supplied spatial graph and with robust finite dimension/scale; G2 quantum dynamics/constraints with independent conserved energy; G3 robust common propagation; G4 same-law backreaction from complete conserved energy/stress; G5 viable tensor response / finite induced stiffness and cutoff. P = within cited model; C = conditional/partial; U = unestablished by selected construction; F = specified counterexample. No aggregate candidate score or generalized theory rejection.

## 2. Discriminant A — a more general limitation than GC32C

Suppose H(w)=H0+sum_e w_e O_e, where the Hilbert space and any constraints are fixed independently of w; all w_e>=0 and O_e are positive semidefinite. This can be an interacting quantum system, including positive link strain, not necessarily Gaussian. Assume stable H, well-defined ground energy and finite partition function where used. Then

    E0(w) = inf_normalized psi <psi|H(w)|psi>
    F_beta(w) = inf_density rho [Tr rho H(w) + beta^(-1) Tr rho log rho].

Each objective for a *fixed* psi or rho is affine in w. An infimum over the same domain is concave. Furthermore increasing a particular w_e gives a nonnegative change of every trial expectation, so E0 and F_beta cannot decrease. At differentiability, partial_e E0=<O_e>_0 >=0 and partial_e F=<O_e>_beta>=0. The inequalities, not differentiability, are the general statements. This also works in a fixed conserved-number sector or with a fixed chemical potential term included in H0. Degenerate ground states may produce nondifferentiability. Need not be strictly increasing.

Consequences:
* The empty-link configuration is a minimizer of either matter potential when all links can be independently decreased. Unlike GC32C's strict oscillator proof, need not be unique if some O_e have a zero-expectation state.
* Finite-temperature quantum matter entropy alone cannot reverse this positive-link monotonicity. Entropy of *different link configurations*, with a graph measure and link Hamiltonian, is a different problem.
* Fixed budget plus concavity often rewards uneven response; it does not prove every discrete graph minimum is clustered (and fixed-degree discrete graphs are already extreme points of their continuous polytope).
* Simply adding other positive semidefinite species or interactions of this form cannot fix the unconstrained selection issue. This extends the prior observation to non-Gaussian matter, conditionally.
* Escapes include changing the physical degrees/constraints as links change, independently motivated signed/indefinite link terms, noncommuting dynamical link operators instead of externally imposed classical w, collective kinetic terms and configurational entropy. These avoid the premise; none ensures extended locality. If a quantum-link H still depends on external parameters multiplying PSD operators, the same theorem applies to those external parameters. An arbitrary negative-energy ghost is not an acceptable stabilization mechanism.

## 3. Discriminant B — a competing stable fermionic contribution

Do not insert a negative bosonic zero-point sector by hand. Instead take a conventional finite, number-conserving spinless fermion hopping Hamiltonian H_f=-t sum_ij A_ij c_i^dag c_j with t>0, at half filling (N even, N/2 particles). Its filled-band ground energy is bounded below at finite N and fixed degree. A is symmetric adjacency. The link addition operator is indefinite (unlike the positive Laplacian bilinear). One could equivalently add 2t per particle on d=2 graphs without changing the comparison. No relativistic Dirac-sea subtraction or EM ZPF equivalence claimed; filling and hopping scale are model inputs.

For two competing graphs with equal vertex count and degree, cycle C_n has adjacency eigenvalues 2 cos(2pi k/n). For C12 vs two disconnected C6, N=12, d=2, six fermions:
    E_f(C12) = -(4+2 sqrt(3))t ~= -7.464102 t
    E_f(2 C6) = -8t.
There is a Fermi-level zero-mode degeneracy on C12, which does not affect its ground energy. Splitting into two C6 gives lower energy. It is NOT a general fermion-clustering theorem: N=8, d=2, half-filled C8 energy=-(2+2sqrt(2))t ~= -4.828427t vs 2 C4 energy=-4t, so C8 wins that comparator. Band shell effects can depend on filling/size, which makes a robust geometry law harder, not easier.

For the tested ordinary boson U_b=(1/2)tr sqrt(mu²I+L), mu>0, d=2 so K=cI-A with c=mu²+2>2. Expand convergently:
    sqrt(cI-A)=sqrt(c)[I - sum_{r>=1} a_r (A/c)^r],  a_r>0.
Trace A^r counts closed walks. At each cycle vertex, a length-r walk consists of left/right steps and closes if the signed displacement is a multiple of n. Every walk closing on C12 closes on C6; additional ones occur at r=6 (winding +/-6). Equal N implies tr A^r(2C6)>=tr A^r(C12), strict at r=6. With all expansion coefficients negative beyond the constant, exactly
    U_b(2C6)<U_b(C12) for every mu>0.
Small arithmetic illustrations (whole graphs, same N/d):
| mu | U_b(C12) | U_b(2C6) | difference |
|---:|---:|---:|---:|
| .1 | 7.6705636084 | 7.5823438783 | .0882197302 |
| 1 | 10.0656586065 | 10.0644951022 | .0011635042 |
| 5 | 31.1662090819 | 31.1662090786 | 3.3239e-9 |

Thus for ANY nonnegative sector multiplicities/weights b,f with b+f>0 and t>0, b U_b+f E_f has lower energy on the disconnected 2C6 than C12. No tuning their relative strengths can rescue this specifically tested extended cycle. This is an exact comparator at all positive regulators, not a global minimizer proof or negative verdict on all boson-fermion phases, other fillings, supersymmetric spectra, noncommuting link dynamics or interactions. It is not an attempt to derive the physical ZPF by subtracting arbitrary spectral signs.

## 4. Discriminant C — actual quantum geometry–matter conversion

Graphity's prior work introduced edge quantum degrees; Hamma et al. additionally made link changes exchange with matter and included hopping. Hamma's full dynamics is not equivalent to the following deliberately stripped minimal control. Three abstract sites, one qubit per unordered pair, hard-core bosons at sites. Let

    H = epsilon N_m + Delta N_e -
        g sum_(i<j)(a_i^dag a_j^dag sigma^-_ij + h.c.)
    Q = N_m + 2 N_e.

For epsilon,Delta,g real, g>0, H is Hermitian and time independent, so energy conserved. The conversion term removes one link and creates two bosons or inverse; [H,Q]=0. Q is a conserved occupation resource (not physical energy, even at bare resonance, because exchange contributes to H). It is not a gravitational charge. No graph geometry is given; all pairs eligible.

Exact Q=2 sector: only one occupied link with no bosons, |L_ij>, or two bosons at distinct sites with no links, |P_ij>. For each of the three pairs an independent two-level block:
    H_pair = [[Delta, -g],[-g, 2epsilon]]
    E_± = (Delta+2epsilon)/2 ± sqrt((Delta-2epsilon)^2/4+g²).
If d=Delta-2epsilon, link probabilities in energy eigenstates are
    p_± = partial_Delta E_± = 1/2 ± d/(4 sqrt(d²/4+g²)).
Starting in one pair block cannot move to another: there is no matrix element between distinct pairs. True link/matter superposition and oscillation, but neither a mechanism for growing connected geometry nor native mobile particle. At exact resonance d=0, both eigenstates have link probability 1/2 and identical <N_m>=1 while their energies differ by 2g. They differ in coherent conversion energy. Link detuning response misses that energy exactly. A coherent conversion Hamiltonian can therefore conserve energy yet its simplest matter/geometry observable can be blind to a part of it. No claim larger resource sectors or full Hamma hopping model obey this exact block decomposition; need additional terms/dynamics to escape.

Other quantum link kinetics: for a connectivity graph of all allowed **graph configurations**, choose H_RK=sum_{G~G'} |G-G'><G-G'| (notation |G-G'>=|G>-|G'>), one term per unordered move. It is positive semidefinite and annihilates the equal-amplitude superposition on each connected move component. If all configurations of a fixed edge budget communicate, unique ground over that component is a graph lottery (uniform labeled probabilities), with no preference for manifolds relative to their multiplicities. If connectivity enforced or graph labels identified as gauge, the counting measure must be stated and is itself input. This is a useful null comparable to earlier equal-rate controls, not a theorem against general quantum graph dynamics. Vertex permutation invariance does not forbid an ordered low-dimensional relational geometry: label-averaged adjacency may be constant even for a superposition over all relabelings of a single lattice. Must examine relational invariant distribution of geometry, not individual labels. Symmetry breaking, correlations, and actual dynamics are separate.

## 5. Discriminant D — source identity and a route toward universality

Let H(x) be differentiable in one homogeneous vacuum geometric variable x within a fixed physical Hilbert subspace P; define E0(x) independently from a differentiable vacuum reference. For stationary eigenstates, Hellmann-Feynman gives variation of excitation energy E=E_state-E0 as <partial_x H>-partial_x E0. Define the integrated excess source operator
    S_x=P [partial_x H - E0'(x) I] P.
If this one geometric mode is claimed to couple proportionally to total excess energy in ALL normalized physical states, the precise condition is
    S_x = kappa P [H-E0 I] P
for one kappa independent of state/composition. Necessity from equality of Hermitian quadratic forms: evaluating arbitrary normalized superpositions fixes both real and imaginary off-diagonal elements; sufficiency by taking expectations. If only stationary energy eigenstates are required, this is weaker: each energy degenerate block of S must be kappa(E-E0) times identity, while off-diagonal between nondegenerate levels untested. General local curved spacetime responses couple to the stress tensor, including momentum flux and pressure, not simply to a scalar E for every deformation. The identity is a targeted uniform lapse/clock or proposed *total-energy* source test; not a required equality for arbitrary spatial shear, curvature or volume change, and not Einstein equations.

Application:
* GC32C link stiffness derivative = (q_i-q_j)^2/2. Uniform massive oscillator quantum has no difference across any bond, yet nonzero excess E. The integrated link response cannot equal kappa>0 times H-E0 in all states.
* Q=2 conversion above, at resonance, partial_Delta H=N_e; both split levels have same derivative even though excitation gap 2g. Fails even stationary test for nonzero kappa. Coherent interaction energy is invisible under that variation.
* Control: if H(x)=alpha(x) H_ref + f(x) I with alpha>0, derivative after subtracting vacuum is exactly (alpha'/alpha)(H-E0). All kinetic, rest, interaction, binding contributions scale together, independent of arbitrary species. This is an imposed uniform lapse/common clock scaling, not deriving x, a tensor mode, attractive force, or locally observable effect if only overall time units change. It identifies the kind of operator-level relation a proposed vacuum mechanism must supply, and the terms missed by partial Hamiltonian couplings.
* General target: start with microscopic action S_m[X] and infer G[X] independently from propagation (GC32B). Metric stress defined by delta S_m= -(1/2) integral sqrt(|G|) T_mu,nu delta G^mu,nu on shell, plus possible other field variations (sign convention defined). Check whether variation with respect to fundamental variables really factors through G for all sectors, and whether accessible delta G spans the physical metric polarizations. A universal cone by itself says only about highest derivatives, not full mass/binding response, nonmetric forces, or Einstein-Hilbert field dynamics. An independent energy ledger and projection/Ward identities are needed.

## 6. Literature screen and new alternative

| Candidate construction | G1 | G2 | G3 | G4 | G5 | Why retain / why not sufficient |
|---|---|---|---|---|---|---|
| Quantum Graphity (Konopka et al. 2008) | C | P | U | U | U | Competing valence/loop energies and quantum graph space; reported ordered low-dimensional phase evidence. Valence/loops and scales are input; tuning desired geometry is not deriving it. PI explored previously. No repeat. |
| Quantum matter-link conversion (Hamma et al. 2010) | C | P | C | C | U | Full hopping + exchange on dynamical quantum graph, finite propagation studied; geometry can react to matter; no demonstrated universal stress/metric sourcing. Our minimal control is narrower and fails source criterion. |
| String-net condensation (Levin/Wen 2005) | U | P | C | U | U | Collective topological order gives emergent gauge/fermionic excitations; locality/lattice input, gauge structure is not a spin-2 metric. Topology alone not same speed (GC30B2). |
| Holographic QEC / entanglement equilibrium (Almheiri/Dong/Harlow, Jacobson) | C | C | C | C | C | Encoding and local entropy dynamics can support geometric consistency within chosen duality/geometric assumptions. Not from arbitrary quantum states. GC32A prevents inferring cone from state alone. |
| Positive matter link potentials and simple competing free bands | F | P | U | F | U | Analytical A/B limitations for specified assumptions and comparator; fixed link classical potentials, no geometry or common stress law. |
| Noncommutative Yang–Mills matrix models, IKKT-type | C | C | C | C | C | Geometry encoded by algebra/noncommutativity; scalar/fermion kinetic principal metrics arise from same matrices; conditional one-loop curvature stiffness. Background choice, Lorentzian signature, stable dimension/scale, physical energy/clock, extra modes and complete sourcing still unresolved. This is the most concrete different mechanism to test next, not a full pass. |

Here C at G2 for IKKT is deliberate: the full model is a reduced matrix action/integral, not a demonstrated Hamiltonian time/energy ledger within this screen. Do not automatically assume matrix dimension equals physical spatial dimension. String-net protects topological phase, not necessarily all low-energy speeds. No candidate passes all gates. Weinberg-Witten premise still matters: ordinary Lorentz-covariant stress tensor of an underlying local flat-space QFT with emergent graviton is restricted; pregeometric matrix variables and gravitational gauge symmetries may escape the premises, but do not prove a healthy tensor sector. See GC32A.

### 6.1 Why matrix variables are attractive: one physical differential structure

Unlike a link lattice with multiple independent hopping rates, start with D Hermitian matrices X^a (matrix size N, target index a carries an input target bilinear eta), where the same inner derivations act on matrix fluctuations. In the familiar noncommutative-brane semiclassical regime, matrix commutators approximate Poisson derivatives:
    [X^a, Phi] ~ i theta^{mu nu} (partial_mu x^a)(partial_nu phi).
Define embedding metric g_mu,nu = (partial_mu x^a) eta_ab (partial_nu x^b) and e^{a nu}=theta^{mu nu} partial_mu x^a. Up to volume/conformal factors, scalar action Tr [X^a,Phi][X_a,Phi] has principal tensor
    gamma^{mu nu}=eta_ab e^{a mu} e^{b nu}
                 =theta^{mu alpha} theta^{nu beta} g_alpha,beta.
For a fermion with Dirac matrix action D=Gamma_a ad(X^a), the first-order principal symbol is Gamma_a e^{a mu} k_mu; Clifford anticommutation gives its square = gamma^{mu nu} k_mu k_nu I (irrelevant i/sign convention). Derivative/connection effects are subprincipal. Hence a shared characteristic cone for these two independently different kinetic operators is an algebraic consequence on the same nondegenerate semiclassical background, not equality by speed fitting. The finite-matrix identity D² = eta_ab ad(X^a)ad(X^b) + (1/2) Gamma_ab [ad(X^a),ad(X^b)] with Gamma_ab=[Gamma_a,Gamma_b]/2 is the underlying relation (indices contracted appropriately). The last commutator is ad([X^a,X^b]) and represents field-dependent terms. At finite matrix size this does not itself define a Lorentzian cone; the differentiable brane limit must exist. When theta is nondegenerate, gamma has the signature of g by matrix congruence. In a degenerate or non-Lorentzian background no physical shared light cone follows.

This is a known result of matrix-model geometry (not a new derivation claimed by Model Gc); the significance is that it attacks GC30/31's repeated two-cone obstacle by a different algebraic structure, with a single X entering both sectors. Effective gauge fields can also share the principal metric in suitable backgrounds after removing constraints. Fermion connections, measure/dilaton, masses from extra commutators, gauge terms and radiative corrections can still cause nonmetric/multisector response. Merely inheriting this common operator isn't proof of full stress coupling.

### 6.2 Sakharov connection and limits
The IKKT model originated as a reduction of ten-dimensional super Yang-Mills matrix degrees. Lorentzian simulations by Kim, Nishimura and Tsuchiya reported expanding three spatial dimensions in one regulated formulation; target index dimension/signature and regulators are inputs, and this is evidence under that formulation, not generic stable 3+1 geometry derived here. A simple positive Euclidean classical commutator-squared potential is minimized by commuting matrices, so one must not assume a noncommutative brane is an unconditional stable vacuum. Fermion/quantum contributions, constraints and dynamics matter; Lorentzian action is not a positive Euclidean energy minimization.

Steinacker's 2023 primary calculation of the IKKT one-loop action on quantum spacetime branes with compact fuzzy extra dimensions includes an Einstein-Hilbert term, with Newton scale determined conditionally by the internal Kaluza-Klein scale. Additional torsion, higher-spin and nonstandard terms appear; a later 2026 paper explicitly studies modified gravity at large distances instead of assuming Einstein dynamics exact. This is a direct structural analogue of Sakharov's proposed vacuum-induced stiffness, but one-loop on particular backgrounds and not an empirical validation, proof of standard GR, prediction of G from first principles, or natural particle cutoff. Need an independent reason for the compact scale and phase stability; extra dimensions and embedding provide assumptions. GC29A does not fix this cutoff. The promising part is the coexistence of algebraically shared kinetic geometry and calculable induced action, whereas GC32C tried to select a supplied link spectrum using positive ground energy alone.

## 7. Recommendation — EX-GC33B (proposal only)
Do not run large IKKT simulation or fit another graph Hamiltonian. **Common metric and complete source from matrix geometry**: first do a short algebraic benchmark on a simple semiclassical noncommutative 3+1 brane (positive control with known background, not emergent selection), with scalar, spinor and gauge kinetic terms drawn from one matrix action. Independently infer principal gamma and metric volume factor from each, then vary the same X (or theta and embedding perturbations generated by X). Check full action variation factors through one stress tensor including mass/coupling/binding terms; separate dilaton, torsion and nonmetric variations. Check rank of accessible delta G modulo coordinate/gauge redundancy against required metric perturbations, and compare actual propagating tensor sector rather than count all independent matrix entries as gravitons. For instance delta gamma=delta theta g theta^T + theta delta g theta^T + theta g delta theta^T; not every metric perturbation need be in that image. Independently track energy using a given reference time on this conditional control; do not count supplied clock as emerged. Include constant Poisson background/gauge-off controls, variations in independent kinetic normalization, and explicit scalar mass or internal KK-energy perturbations, not just massless characteristics.

Pass gate is a demonstrated common principal metric AND full energy/stress variation through that metric on the selected physically allowed sector; separately mark restricted metric image, extra forces and energy gaps. Existing one-loop action helps select where to test but is not a pass. Only if the coupling is robust should a phase-selection or induced-stiffness calculation follow; ask whether this framework supplies a mechanism better than GC30/31. If coupling is nonuniversal or lacks necessary metric polarizations, stop that arm, don't optimize the model to GR. This is a proposed new branch; no production authorization assumed.

## Primary literature
- Konopka, Markopoulou & Severini, [Quantum Graphity: a model of emergent locality](https://arxiv.org/abs/0801.0861), Phys Rev D 77 104029 (2008).
- Hamma, Markopoulou, Lloyd, Caravelli, Severini & Markström, [A quantum Bose-Hubbard model with evolving graph as toy model for emergent spacetime](https://arxiv.org/abs/0911.5075), Phys Rev D 81 104032 (2010).
- Levin & Wen, [String-net condensation: a physical mechanism for topological phases](https://arxiv.org/abs/cond-mat/0404617), Phys Rev B 71 045110 (2005).
- Almheiri, Dong & Harlow, [Bulk locality and quantum error correction in AdS/CFT](https://arxiv.org/abs/1411.7041), JHEP 04 163 (2015).
- Jacobson, [Entanglement equilibrium and the Einstein equation](https://arxiv.org/abs/1505.04753), PRL 116 201101 (2016).
- Ishibashi, Kawai, Kitazawa & Tsuchiya, [A Large-N Reduced Model as Superstring](https://arxiv.org/abs/hep-th/9612115), Nucl Phys B 498 (1997).
- Kim, Nishimura & Tsuchiya, [Expanding (3+1)-dimensional universe from a Lorentzian matrix model](https://arxiv.org/abs/1108.1540), PRL 108 011601 (2012).
- Steinacker, [Emergent Gravity from Noncommutative Gauge Theory](https://arxiv.org/abs/0708.2426) (2007); [Emergent Geometry and Gravity from Matrix Models](https://arxiv.org/abs/1003.4134) (2010, review by researcher for principal action details); [One-loop effective action and emergent gravity on quantum spaces in the IKKT matrix model](https://arxiv.org/abs/2303.08012), JHEP 05 129 (2023); [Modified gravity at large scales on quantum spacetime in the IKKT model](https://arxiv.org/abs/2601.08031), JHEP 04 044 (2026).

**Disposition:** paper screen complete. General positive-link theorem and exact competing-band/conversion counterexamples, source-identity benchmark and conditional alternative identified. No physical gravity mechanism established. Review/PI closure pending; GC33B proposed only.
