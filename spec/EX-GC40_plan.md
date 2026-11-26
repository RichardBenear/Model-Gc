# EX-GC40 — geometry as counts: entanglement equilibrium and event counting on a pair network

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread ("yes, go ahead and write the plan", 2026-11-26 15:18 UTC). Written before any run. **Running needs a separate go-ahead.** Proposed cap 4 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why this study

**Source of the idea.** The Chain Theory compendium (working draft, Oct 2026, supplied by the PI) starts from one axiom: nothing exists but finitely many maximally entangled pairs of finite-dimensional systems.
- We take only that axiom.
- We do not use its link genres, its webs, its constants, or its gravity results. Its own field equation fails at all 60 cells it tests (its Open result 92.1), so there is nothing there to reproduce.

**What our own results say.**
- **39A:** on a lattice with weights on fixed links, entanglement equilibrium does not single out a geometry. Both responses are boundary (cutoff) effects that scale as ℓ². The answer also depends on how "same volume" is defined.
  - The proposed fix was that **area should be the number of cut links and volume the number of sites**, which is exactly a pair network.
- **37B:** clocks are one-site events (rate N) and signals are two-site coincidences (rate N²). Together they give the weak-field metric with bending ratio 2. This was shown only on a regular lattice, where stencil dispersion added +9.6% to the absolute bending.

**The two questions.**
- **40A:** on a pair network, does entanglement equilibrium give a Newton-like law, with a fixed sign and a G that is independent of block size?
- **40B:** does 37B's counting rule still give an isotropic light speed and GR's factor-2 bending on an irregular network, where there is no lattice direction?

## 40A — entanglement equilibrium on a pair network

### Model

**Substrate (a decorated cubic lattice of free fermion modes).** It is Gaussian, so 3D sizes are cheap.
- Periodic cubic lattice, L sites per side.
- Each site v carries one fermion mode per link end: six modes a_{v,±x}, a_{v,±y}, a_{v,±z}.
- **Pair term:** H₀ = −Σ_links (a_{v,+a}† a_{v+a,−a} + h.c.). Each link joins its own two modes, one at each end.
- At half filling, the ground state of H₀ is a product of one bonding orbital per link. **Each link is then one maximally entangled pair (one ebit),** which is the axiom's vacuum.

**Dressing (moving away from pure pairs).**
- H = H₀ + λH₁, with H₁ = Σ_v Σ_{l<l'} (a_{v,l}† a_{v,l'} + h.c.), an on-site coupling among the six modes of each site.
- λ = 0 is the pure pair network. As λ grows, entanglement spreads beyond single links.
- Values: λ ∈ {0, 0.1, 0.2, 0.3, 0.5}.
- The band gap is computed for each λ and recorded. Gapless values are reported but kept out of P3.

**Larger link registers.**
- A variant with k = 2 modes per link end (link dimension d = 4) is used for the area check only.
- A qutrit (d = 3) is not natural for fermion modes, so it is not included.

**Blocks.**
- Cubes of ℓ = 2–6 sites, all modes included.
- Volume V = ℓ³ sites. Area A = number of cut links = 6ℓ².

**Numerics.**
- Correlation matrix C = ⟨a†a⟩ from the 6-band Bloch Hamiltonian by FFT (L = 16).
- S_B from the eigenvalues of C restricted to B: S = −Σ[c ln c + (1−c) ln(1−c)].
- Real-space dense diagonalisation (L = 10, 6,000 modes) for the rewiring test, which breaks translation invariance.

### Readouts

**1. Area coefficient.** s_link(λ) = S_B / (cut links), for each ℓ. Expect ln 2 at λ = 0 (ln 4 for d = 4).

**2. Geometry response (rewiring at fixed volume).**
- Take two cut links (a,b) and (c,d), with a and c inside B and b and d outside. Rewire them to (a,c) and (b,d).
- This removes 2 cut links. Every site keeps its degree and its modes, so volume is exactly fixed.
- η_count = −δS_B/2, compared with s_link.
- Positions tested: face centre, near an edge, and next to a corner; ℓ = 3, 4, 5.
- In a count geometry this is what "changing the area at fixed volume" means. The weight rescalings of 39A have no counterpart here.

**3. Matter response (k_B).**
- A first-order Gaussian deformation of the vacuum, δC = iε[F, C], the fermionic analogue of 39A's squeezing.
- F mixes the highest filled and lowest empty band over long wavelengths (|k| ≤ k₀ = 2π·2/16), paired with k + q, q = 2π/16 ẑ.
- The local energy density δρ(x) is measured from the deformed state, not assumed.
- k_B = δS_B / δρ(c), with δρ the block average, as in 39A.
- **Scaling with ℓ.** Jacobson's bulk first law needs k_B ∝ ℓ⁴ (the boost-weighted energy of a ball). 39A found ℓ². Fitted exponent over ℓ = 3–6.

**4. Equilibrium G.**
- δS_geo + δS_matter = 0, with δS_geo = s_link × δ(cut links). This is 39A's chain with η_eff replaced by s_link.
- With 39A's convention, G_eq(ℓ) = 3k_B / (4π s_link ℓ⁴).

### Predictions

- **P1 (area law, structural):** at λ = 0, S_B = ln 2 × (cut links) to 1e-10 for every block, and ln 4 × (cut links) for d = 4. For λ > 0, s_link decreases with λ.
- **P2 (geometry is counts):** η_count / s_link = 1 within 10% for every λ, ℓ and rewiring position.
  - At λ = 0 this holds exactly.
  - For λ > 0 it is the real test: is the entropy carried by each cut link the same everywhere?
- **P3 (pure pairs do not gravitate):** at λ = 0, δS_B has no first-order part: |δS_B|/ε falls at least linearly as ε goes from 1e-3 to 1e-4, so k_B = 0.
  - Reason: inside B every link is a pure orbital (c = 0 or 1, whose eigenvalues move only at second order), and every cut mode sits at c = 1/2, the maximum of the entropy. So no first-order change of the state moves S_B.
  - Consequence: a vacuum made only of maximally entangled pairs has no first law. Matter cannot act on geometry through entanglement there.
- **P4 (λ > 0, the main test):** G_eq(ℓ) has one sign for all ℓ and λ, and is independent of ℓ within 20% over ℓ = 3–6 (gapped λ only).
  - Expectation: uncertain. If k_B scales as ℓ², as in 39A, then G_eq ∝ ℓ⁻² and P4 fails.
- **P5 (descriptive):** the sign of G_eq (attractive or repulsive), and G_eq against 1/(4 s_link). In Jacobson's continuum derivation these two are equal.

### Checks

- **G1:** FFT correlators agree with dense real-space diagonalisation at L = 6 to 1e-12. λ = 0 gives exactly one ebit per link.
- **G2:** linearity, for λ > 0. δS_B from ±ε agree in magnitude to 1%. (At λ = 0 the response is second order by P3, so G2 does not apply.)
- **G3, separate path for k_B:** the first law δS_B = Tr(K_B δC_B), with K_B = ln((1 − C_B)/C_B) computed in arbitrary precision (python-flint, as in 39A). It must equal the finite-difference δS_B to 1% for ℓ = 2 and 3.
- **G4:** the rewiring δS from the dense solve is reproduced by a second route: a rank-4 update of the ground state from the λ-perturbed Green's function. Agreement to 1e-8.

## 40B — event counting on an irregular pair network

### Model

**Graph.**
- 2D Poisson sprinkling at density 1, in a disk of radius R = 150 (about 70,000 sites).
- Links join every pair of sites closer than r_c = 1.5 (mean degree about 7). Each link is one pair.
- No lattice directions exist, so there is no stencil dispersion.

**Firing and signals.**
- Each site fires as a Poisson process at rate N(x).
- **Two-site rule (37B):** a signal crosses link (i, j) when one end fires and the other fired within τ = 0.01.
- **One-site control:** a signal crosses when the sending site fires.

**Arrival times, two versions.**
- **Deterministic:** each link's time is its mean crossing time (the 37B rate sum). Arrival time T(x) is the weighted shortest path (Dijkstra).
- **Stochastic:** each link's time is one sample of the 37B `walk` rule. T(x) is the first-passage time, over 20 graph-and-firing realisations.

**Lapse.**
- Prescribed: N = 1 − 2GM/√(r² + r_s²), with GM = 2 and r_s = 3.
- Where N comes from is not studied here. 39C showed that it must not be a free radiating field.

### Readouts and predictions

- **P6 (isotropy):** with flat N = 1, the speed d/T at d = 100, averaged over realisations, agrees across 16 angular sectors to 1% (deterministic and stochastic).
- **P7 (speed exponent):** with uniform N from 0.5 to 1, v ∝ N^p, with p = 2.00 ± 0.05 (two-site) and 1.00 ± 0.05 (one-site), for both versions.
- **P8 (Shapiro delay ratio):** source and receiver at x = ∓120, with impact parameters b = 20, 40, 60. The ratio of the two-site to the one-site delay tends to 2, within 3% at b = 60 (deterministic). The stochastic ratio is reported, descriptive.
- **P9 (bending):**
  - The deflection is the turning angle of the shortest path from the source to far receivers, from the Dijkstra predecessor tree.
  - Two-site / one-site ratio: 2 within 5% at b = 40 and 60.
  - Absolute two-site deflection: within 5% of the continuum eikonal for N⁻², after the flat path-stretch factor is calibrated out (deterministic).
  - This is the test that 37B's +9.6% excess was stencil dispersion.

### Checks

- **G5:** on a fine square lattice, the deterministic Dijkstra reproduces 37B's rate-sum delays to 1%.
- **G6:** the mean stochastic link time against the rate-sum formula, including 37B's known 0.84% start-up offset.
- **G7:** two graph seeds give P8 and P9 values that agree within their stated tolerances.

## 40C — follow-on, after 40A: is a shared-entanglement "hub" a mass? (outline only)

The compendium builds mass from reads that pass through hubs, systems shared by more than one link (its Part X, Theorems 185.6–185.9). 40A's uniform λ spreads sharing everywhere, which makes it a vacuum property. 40C puts the sharing in one place.
- **Body:** λ ≠ 0 only on a small cluster of sites (radius 1–2) in an otherwise pure pair network, or on a dressed background at λ₀ with λ₀ + Δλ on the cluster.
- **Readouts:**
  - the energy the hub carries;
  - whether it is a localised, gapped excitation (a mass);
  - the profile of δS_B and of the equilibrium rewiring around it with distance;
  - whether that profile falls off as 1/r (Newton) or stops at the cluster edge.
- 40C gets its own predictions, written after 40A's results and before any 40C run.

## What a result would mean

**40A, if P3 holds and P4 passes at λ > 0.**
- Pure maximally entangled pairs cannot be the whole story, because a vacuum of pure pairs has no first law. Pairs that share entanglement with their neighbours can be.
- Matter would then set the geometry through counts, and the rewiring that equilibrium asks for would be the field equation.

**40A, if P4 fails as in 39A.**
- Counting fixes the geometry side (P2), but the matter side is still a boundary effect.
- Then entanglement equilibrium is not the mechanism on any substrate of this kind. That would close the H-A-through-entanglement line of 38 and 39.

**40B, if P6–P9 pass.**
- 37B's counting rule does not need a lattice: an irregular pair network carries GR's light propagation with no metric put in.
- This would be the strongest form of the event-counting result so far.

## Not in scope

- Chain Theory's link genres, constants and Standard Model claims.
- Dynamics of the graph itself: no rewiring in time, no source law for N.

## Outputs

`src/exgc40a.py`, `src/exgc40b.py`, `results/EX-GC40A/`, `results/EX-GC40B/`.
