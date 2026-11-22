# EX-GC31A — Does causal invariance enforce a common metric?
2026-11-22 | ChatGPT | EXPLORATORY | PI authorized | Analytical test, no production simulation; negligible compute relative to 2 core-hour cap. Review/PI closure pending.

Registered design: [plan](../../spec/EX-GC31A_plan.md). GC30A/B/B2/B3 remain closed. The question is Claude's suggestion: can a pregeometric causally invariant rewriting system force all excitations onto one effective light cone? This tests a logical implication, not a particular unrestricted ModelG rule.

## 1. Operational gates
H1: same event dependency DAG under all fair schedules of a specified bounded system. H2: common structural reachability bound. H3: same sector-resolved exact retarded support/front. H4: common longwave characteristic cone in one independently reconstructed clock/space calibration. H5: manifoldlike Lorentzian order and volume/clock normalization. H6: dynamical energy sourcing and universal massless tensor coupling. P = established within model, F = explicit failure, U = unestablished. H4 does not require equal velocities for massive particles.

## 2. Exact counterexample A: one rewrite, two ballistic speeds
Construct a single uniform deterministic port rewrite, on an infinite line or ring (before wraparound), with finite depth T. Event (t,x) consumes one-use inputs from left, right and self predecessor events and outputs unique tokens on its left, right and self successor ports. Initial tokens at layer zero are supplied; each event fires once when all ports are ready. The labels index the given combinatorial wiring and are not scheduler time or an assumed spacetime metric. Every ready event has disjoint input tokens from every other ready event. Induct on DAG depth: each event receives uniquely determined inputs from its predecessors, produces uniquely determined outputs and dependency edges, and is enabled once they exist. Every fair complete execution therefore differs only by a linear extension of the same labelled DAG; every partial schedule extends to the same complete cut. This establishes exact scheduling invariance, not merely agreement of species counts. There is no matching ambiguity in this restricted port system; it does NOT prove causal invariance for unrestricted unlabelled hypergraph rewrites with competing overlaps. The wiring is given, not dynamically generated.

Each message carries the appropriate channels of F and S. Use the single local rule (right-moving sector shown; left obtained by x→−x):

    F(t+1,x) = F(t,x−1)
    S0(t+1,x) = S1(t,x−1)
    S1(t+1,x) = S0(t,x)

The F channel and the S1→S0 channel use the right port, S0→S1 uses self; the reflected channels use left and self. Thus every type of edge has a real dynamical dependency; no padding/dummy edges. The global update on all channels is a reversible permutation, conserves their total content, and is translation and reflection symmetric. An initial point signal travels t links for F, floor(t/2) for S0, ceil(t/2) for S1. The S signal advances one link in two layers; on the two-step stroboscopic description its gapless frequency slope is 1/2, versus 1 for F (common dimensionless units). The maximal graph reach is ±t, but operational sector fronts are 1 and 1/2. Neither a different coordinate calibration nor a mass gap accounts for the difference.

H1 and H2 pass; H3 and H4 fail. Even exact deterministic scheduling independence and one rule do not enforce universal sector propagation.

## 3. Exact counterexample B: common microscopic front and exchange symmetry still insufficient
Use the same deterministic event architecture, transporting current and previous values as required (history on the self port), with one two-component local wave rule:

    U(t+1,x) = 2 U(t,x) − U(t−1,x)
              + K [U(t,x+1) − 2U(t,x) + U(t,x−1)]
    K = [[1/2, 1/4], [1/4, 1/2]]

This is invariant under swapping the two components. The symmetric and antisymmetric eigenchannels u± = (U1 ± U2)/sqrt(2) have eigenvalues λ+ = 3/4, λ− = 1/4. Fourier substitution gives exactly

    sin²(ω±/2) = λ± sin²(k/2).

All nonzero k have real stable frequencies for these λ in (0,1); k=0 has the usual constant and linear-in-time zero-frequency wave solutions. Both branches are gapless. Their longwave speeds in identical units are sqrt(3)/2 = 0.866025 and 1/2 = 0.5 (ratio sqrt(3)); the principal spatial coefficients differ by a factor of 3.

For an impulse in either eigenchannel at t=0, with u(−1,x)=0, no response is possible at |x|>t. At x=±t there is exactly one contributing outward path, hence u(t,±t)=λ^t (unit initial impulse). Both exact microscopic support fronts are therefore speed 1, with nonzero extremal amplitude at every finite depth before wraparound.

| t | symmetric front amplitude | antisymmetric front amplitude |
|---|---:|---:|
| 10 | 0.0563135 | 9.53674e−7 |
| 40 | 1.00566e−5 | 8.27181e−25 |

These are exact powers evaluated for illustration, not fitted or thresholded front velocities. A shared maximum microscopic reach and even shared exact support are NOT the same as a shared effective characteristic speed. With lattice spacing a, layer interval Δt, and fixed a/Δt, taking the longwave continuum limit gives ∂²_t u± = (a/Δt)² λ± ∂²_x u±. Extremal microscopic amplitudes vanish exponentially with layer count in that limit. The two resulting wave cones are distinct. Exchange symmetry allows the off-diagonal kinetic term; it does not force equal eigenvalues.

H1, H2 and H3 pass (H3 here means the common extremal front); H4 fails. H5/H6 are not demonstrated in either counterexample. No claim that the prescribed wiring emerged from the rule.

### Positive control and criterion
Replace K by λ I with λ=1/2, and allow distinct diagonal lower-order mass terms −(m_i Δt)² U_i. The exact dispersion is sin²(ω_i/2)=λ sin²(k/2)+(m_i Δt)²/4, stable in the range λ+(m_i Δt)²/4<1. In the continuum, Ω_i² = c² p² + m_i², c=(a/Δt)sqrt(λ). Distinct massive group velocities are compatible with a shared characteristic cone. Both microscopic extremal fronts remain speed 1. A common principal symbol passes H4 by construction, not emergence. A useful sufficient condition for unconstrained propagating fields is a factorization P_ab(p)=Z_ab g^{μν}p_μ p_ν, with nonsingular kinetic normalization Z and one common metric; constrained systems require projection to physical modes. The research challenge is derivation of this property from native rules, not assuming it in the update.

| Arm | H1 | H2 | H3 | H4 | H5 | H6 |
|---|---|---|---|---|---|---|
| A: delayed channel | P | P | F | F | U | U |
| B: exchange-symmetric wave | P | P | P (front) | F | U | U |
| Control: common kinetic term | P | P | P (front) | P (imposed) | U | U |

## 4. Candidate screen
* Causal sets: a manifoldlike causally distinguishing Lorentzian order constrains the conformal metric; volume/counting can supply scale under appropriate sampling assumptions. Not every poset is such an order. Even given that geometry, species' retarded operators are additional dynamical information. Johnston's construction of causal-set propagators is a concrete example of needing a propagation prescription. Promising geometry reconstruction, not an automatic universality theorem.
* Native causal-invariant rewriting: best alignment with original ModelG. Exact equivalence of causal histories, not just isomorphism of end-state counts, is the essential first gate. Neither counterexample excludes a more restrictive native rule that derives a unique kinetic operator and whose topology grows dynamically. Need demonstrate that stronger property, including competing critical overlaps. The examples here are deliberately weaker port-rewriting constructions, not an assertion about all Wolfram-style rules.
* Quantum causal networks and entanglement: a shared quantum causal circuit can bound information transport; the excitation spectrum depends on its dynamics. Static entanglement alone cannot determine operational rates: take independent local Hamiltonians H1+H2 and H1+b H2, b>0. For unique sector ground states the product ground state is identical, as are its equal-time correlations/entanglement across spatial partitions, while excitation frequencies and propagation rates in sector 2 scale by b in a common clock. This is a logical counterexample to static state data alone fixing every sector's speed, not a refutation of holographic/CFT reconstructions with extra dynamical assumptions. Entanglement is not itself an extra signaling channel. Quantum causal automata provide useful models for deriving, rather than assuming, longwave kinematics.
* Strong interacting fixed point: reserve route from GC30B3; weak attraction to a common cone is not enough over finite scale. A preferred-frame UV sector can regenerate Lorentz violation through interactions (Collins et al., particular perturbative example), so any candidate must address stability of its universality. No general no-go for strong fixed points is claimed.

Parent ModelG check (read-only): [MG189](https://github.com/RichardBenear/ModelG/blob/main/results/MG189_causal_invariance/MG189_results.md) found 1192 disjoint event pairs commute but count reconvergence was not enough: R1 creates G and R2 consumes new G versus initial G, giving identical counts but different causal edges. Full causal completion unresolved; the finest scheduler comparison was statistically inconclusive. Do not promote that result to a particle-scale causal invariance proof. The [Claude GC30 review](../../review/EX-GC30A-B3_claude_review.md) is consistent with this distinction. GC30 stays closed.

## 5. Implications for gravity and proposed next gate
Causal invariance is a useful necessary design target for this proposed route: no physical dependence on arbitrary scheduling. It does not alone imply operational universality, much less an Einstein metric. A structural order potentially supplies kinematic geometry; an independent dynamical law must explain how the relation changes with energy, and why all excitations (including composite matter and gravitational perturbations) respond to it. Sakharov vacuum-induced curvature stiffness remains downstream of deriving a common metric. A genuine pregeometry could evade premises of Weinberg–Witten (e.g. no microscopic Lorentz-covariant local conserved stress tensor); that is not evidence it has a universally coupled massless helicity-2 mode. No gravitational force, Newton constant or equivalence-principle precision is claimed.

Recommend EX-GC31B, only after PI authorization: bounded native-rule qualification, not a large simulation or re-run of MG189. Select a candidate with genuine alternative overlapping matches and no imposed lattice/global scheduler clock. Pre-register:
1. Compare causal DAGs from producer/consumer/read dependencies for exact bounded alternative fair completions or critical pairs (not counts); distinguish a bounded pass from a global theorem.
2. Reconstruct dimension, causal order and interval-volume/clock calibration from the unperturbed structure, independently of particle velocities. Verify if a meaningful continuum scaling regime exists. A small 1+1 example is a diagnostic, not gravity.
3. Identify at least two genuine propagating, preferably gapless perturbations of the same native rule (including background and defect perturbations if physically available). Compare paired intervention histories on relationally aligned events; if graphs change and alignment is ambiguous, report unresolved rather than infer a speed.
4. Measure structural reach, exact retarded front, amplitude and dispersion/characteristics separately, across types and scale, in one calibration, without type-specific adjustments. Separate massive group speeds from limiting characteristics; test robustness to topology-preserving/background perturbations and directional anisotropy.
5. Stop if causal completion fails, there are not two propagating modes, or topology/clock is imposed to manufacture equality. A pass only qualifies a candidate for subsequent tests of energy-like backreaction and universal response, with gravitational sourcing and conserved energy/momentum still unestablished.

The essential design change is to search for a rule in which the causal relation controls the dynamics of every physical mode, not merely one on which every mode is allowed to live.

## Primary references
- Bombelli, Lee, Meyer, Sorkin (1987), Spacetime as a causal set, https://doi.org/10.1103/PhysRevLett.59.521.
- Johnston (2010), Quantum fields on causal sets, https://arxiv.org/abs/1010.5514.
- Gorard (2020), Some relativistic and gravitational properties of the Wolfram model, https://arxiv.org/abs/2004.14810; causal set connections https://arxiv.org/abs/2011.12174.
- Bisio et al., Quantum cellular automaton theory of light, https://arxiv.org/abs/1601.04842.
- Van Raamsdonk (2010), Building up spacetime with quantum entanglement, https://arxiv.org/abs/1005.3035.
- Jacobson (2015), Entanglement equilibrium and the Einstein equation, https://arxiv.org/abs/1505.04753.
- Collins et al. (2004), Lorentz invariance and quantum gravity: an additional fine-tuning problem?, https://arxiv.org/abs/gr-qc/0403053.
- Weinberg & Witten (1980), Limits on massless particles, Phys. Lett. B 96, 59–62.

Disposition: analytical tests complete, generic implication H1+H2 ⇒ H4 rejected; native pregeometric hypothesis remains open. No next production authorized.
