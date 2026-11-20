# Gravity, spin and inertia from reaction–diffusion particles: findings of the Model Gc programme to date

**Draft v0.1 for review, 2026-11-20.** PI: Richard Benear. Drafted by Claude (Anthropic) at the PI's request, from the
records in this repository; to be reviewed by ChatGPT/Codex and revised in iterations. Every result below except MG195
is EXPLORATORY under `RESEARCH_PROTOCOL.md` and supports no claim on its own. Model Gc, its variants, the lattice rules
and the couplings are our construction, not Paul LaViolette's. Each statement points to the record it rests on; the
study index in Appendix A lists them all.

**Review notes for the reader.** Sections marked *(interpretation)* go beyond what a record measured. Numbers are
copied from the notes cited; a reviewer should check them against those notes and the raw JSON records.

---

## Abstract

We asked whether particles built as localized reaction–diffusion structures, in the spirit of LaViolette's Model G,
can show gravity, the equivalence principle, spin-½ behaviour and Newtonian inertia, and whether all of this survives
on a discrete lattice of tokens as a step toward Wolfram-style rewriting rules. In Model G itself, a conservation
identity forces a stationary particle's net G disturbance to vanish, so it has no long-range gravity. We therefore built
**Model Gc**, in which a carrier G is conserved and consumed only inside particle cores. This gives a true 1/r field,
and particles that drift toward lower G attract as 1/d²; a registered study (MG195) supported the inverse-square law
at 1–2%. Free fall in Model Gc is not universal across particle types or for bound composites. A recoil rule (each
particle is pushed by the tokens it absorbs, at a fixed rate) restores universal free fall in the weak regime. For
spin, a degree-1 knot of lattice frames under a local three-node rule is stable, keeps its topological charge exactly,
needs a 720° turn to return to itself, and shows the classical fermion exchange signature (Finkelstein–Rubinstein);
the knot can carry a Gc particle that sources the 1/r field and falls at the correct rate. Inertia is where the
programme stalls. Absorbing tokens from a medium at rest produces drag, not inertia, so particles have a universal
terminal speed rather than Newtonian motion. A momentum (flow) field lets the particle move and coast as one body, but
its inertia is assigned, and the inertia the medium itself generates is 20 to 400 times too small against the
absorption drag and scales with the wrong power of the particle's strength. We argue that this is structural: gravity
from net consumption of a conserved carrier requires a medium rest frame, and that rest frame causes the drag. A short
analysis and a 1+1 D check (EX-GC25) show that gravity can survive in a medium with no rest frame if G becomes a field
coupled to particle energy rather than a consumed carrier. The equivalence principle then follows from von Laue's
theorem and inertia comes for free. This would be a change of model, which we put forward for review.

---

## 1. Background and motivation

**Model G and the hypergraph programme.** LaViolette's Subquantum Kinetics describes particles as Turing-type
dissipative structures in a reaction–diffusion "ether" (Model G). The parent repository (RichardBenear/ModelG) mapped
Model G onto discrete token and hypergraph rules and tested, under pre-registered gates, whether the discrete rules
reproduce Model G's particles. The long-term aim is a link to quantum mechanics and gravity, which is why causal
invariance of the rule set matters (Wolfram 2020).

**Why Model G has no long-range gravity.** LaViolette's gravity is a G-well around each particle. In the parent
programme we found that the Model G rate equations imply an identity, (q − g·p)∫φ_G = 0, so a stationary particle has
zero net G disturbance (`docs/lineage/`, EX002–EX004; ModelG memory of the literature check). Growth and variants
(Model G-2, hybrids; EX004–EX007) gave only exponentially screened wells, sometimes repulsive. LaViolette's 1/r gravity
assumes a net exterior G deficit, which the identity rules out for stationary Model G particles.

**The design that followed.** A 1/r field in 3D needs a conserved carrier with no bulk source or loss and a net sink
inside the particle. That is Model Gc (Section 2). Everything after Section 2 is built on it.

---

## 2. Model Gc

Definition: `spec/MODEL_Gc.md`. Gc1 (gravity only):

```
X_t = D_X ∇²X + A − (1+B)X + X²Y + c_x k G W
Y_t = D_Y ∇²Y + B X − X²Y + c_y k G W
W_t = D_W ∇²W − W + c Y W² − W³
G_t = D_G ∇²G − k G W
```

X and Y form a Brusselator ether (A = 140/9, B = 29, D_Y = 12). W is a core species that is zero in the ether
(c·Y₀ = 1.86 < 2) and is sustained inside a particle. G is conserved in the bulk, held at G∞ = 1 far away, and consumed
only where W is present. Consumed G feeds X or Y. Outside the core ∇²G = 0, so G = G∞ − M/(4πr), where M = ∫kGW dV is
the consumption rate. Gc2 adds a conserved signed carrier C and two core types that make or consume it (a charge
analogue; EX010).

---

## 3. Gravity: the most secure result

| Study | Result | Record |
|---|---|---|
| EX008 | Stable 3D particle. Far field G∞ − M/(4πr) to 1e-6, independent of domain size. Linear response: the particle drifts toward lower G, implying 1/d² attraction. | `results/EX008_longrange/` |
| EX009 | Direct axisymmetric pair solve: 0.984–0.986 of the 1/d² prediction at d = 20–30. Short-range Turing repulsion near d ≈ 8. Force ∝ k², 1/r range down to k = 0.001. | `results/EX009_followups/` |
| EX010 | Gc2: tuned mirror charges repel/attract as like/unlike. With G added, every type drifts toward lower G; mobilities within 1–2%. Reciprocity is tuned at one point, not built in. | `results/EX010_charge_gravity/` |
| **MG195** (confirmatory, PI closed) | Inverse-square law registered in advance. P1, P2 SUPPORTED (v/prediction 0.98–0.99, fitted exponent 1.98). C2 (Yukawa control) PASS. P3 UNRESOLVED: a short-range Turing tail set the null floor too high. | `results/MG195_inverse_square/` |
| EX-GC09 | On a 3D lattice (h = 0.5) as whole tokens, the far G deficit matches the lattice point-sink field within 1–4%. The particle survives at Ω = 1000 (M −6%) and 300 (−20%) over T = 30, fades at 100 and dies at 30: about 10⁴ core tokens are needed. | `results/EX-GC09/` |

*(Interpretation.)* Model Gc's gravity is a sink-flow (Le Sage-like) mechanism: the field is a steady inflow of a
conserved carrier toward consumers. This is what gives true 1/r range, and, as Section 6 shows, also what causes the
inertia problem.

---

## 4. The equivalence principle

### 4.1 In Gc1, free fall depends on particle type and binding

| Study | Result | Record |
|---|---|---|
| EX-GC01 | v_i = μ_i M_j/(4πd²) within 1–2%. Action and reaction balance (M_a v_a = M_b v_b) only when mobilities match; a particle with doubled k falls 1.8× faster. | `results/EX-GC01/` |
| EX-GC02 | A slow core species gives real inertia (drag linear in the response time; effective mass confirmed by a direct time response), but the threshold depends on type (17.5–40). | `results/EX-GC02/` |
| EX-GC03 | Active and passive gravitational mass agree within 5%, but drag is set by the shared X/Y core, so fall speeds spread 35%. Tuning a slow G brings this to 4–5%, only by tuning. | `results/EX-GC03/` |
| EX-GC04 | Identical particles bind at d = 7.145 (and 10.11); the bound pair falls at 0.939 (0.972) of a single particle's speed because consumption depends on local G (dM/dG = 4.58). | `results/EX-GC04/` |

### 4.2 The recoil rule restores universal free fall (weak regime)

EX-GC16 replaced G-dependent consumption with fixed-rate absorbers that are pushed along the arrival direction of the
tokens they absorb (a token-lattice simulation, Ω = 1000 tokens per site).
- Five particle types fell at m₀E/3 within 0.2–0.5% (m₀ = 50, local depletion about 1%).
- Two absorbers attracted with equal and opposite forces, matching the lattice mean-field prediction within 1–3% at
  d = 2–4 and 9–11% (1σ) at d = 6–8; the prediction tends to m₀²/(12πΩd²). Active and passive mass are the same number.
- At strong consumption (m₀ = 500, 13% depletion) types fell 2–7% faster than m₀E/3, ordered by how concentrated the
  consumption is, as the mean field predicts. This is a self-energy (Nordtvedt-type) violation that scales as m₀/Ω.

Record: `results/EX-GC16/`.

---

## 5. Spin and statistics

### 5.1 Route 1, a rotating particle: did not work
No non-round core mode of Gc1 oscillates (EX-GC05). A slow core inhibitor ("Gc3") gives a breathing clock, not
rotation (EX-GC06), and in 3D time runs the breathing grows until the core collapses (EX-GC07). Each fix adds another
tuned ingredient. Records: `results/EX-GC05/` to `EX-GC07/`.

### 5.2 Route 3, structure in the lattice: works in principle
- **Belt trick (EX-GC08).** Each node carries one of the 24 rotations of a cube plus an SU(2) "lift" (48 elements).
  A local three-node no-cut rule (short lifts must compose around every triangle) gives a conserved Z₂ index: one turn
  sticks, two turns untangle. Pairwise rules leak or freeze.
- **Tokens on that lattice (EX-GC09)** and **a particle tethered by its own G supply (EX-GC10):** a one-turn twist
  collects next to the core and lowers M by 7.9%; two turns untangle (M within 0.6%). The twist is imposed by hand and
  the coupling strength is chosen.
- **No exchange sign from a tether (EX-GC11).** The frame-field tether has no identity (paths through it can be redrawn
  around a core), so swapping particles changes nothing. A field version of spin–statistics needs the particle itself to
  be a topological soliton.

### 5.3 The particle as a knot (soliton route)
- **A stable lattice knot (EX-GC12).** A degree-1 texture (a "hedgehog" rounded to the 48 lifts) relaxes under the
  three-node rule to a compact knot of about 150 mismatched edges. Its degree, the conserved particle number, stays
  exactly 1. Under a pairwise rule it unwinds. Two same-sign knots keep degree 2.
- **Exchange signature (EX-GC13).** Computed loop classes: 360° turn = 1, 720° = 0, swap of two knots = 1, swap plus
  360° = 0, degree-2 knot turned 360° = 0. All 24 computations agree with Finkelstein–Rubinstein (1968), and each loop
  can be carried out on the lattice by allowed moves. This is the classical (topological) form of the fermion sign.
- **The knot as a Gc particle (EX-GC14, EX-GC15).** W autocatalysis is gated by the smoothed knot density, and the knot
  is pulled toward W (energy −λΣWρ). The knot carries a stable Gc particle that sources the 1/r field (within 0.1% of a
  lattice point sink); without a knot the particle dies. At h = 0.25 the knot particle falls at −4.998×10⁻⁴, identical
  to a plain Gc particle to four digits and within 2% of the continuum value. Without the pull (λ = 0) the knot shrinks
  and the particle dies.

**How the knot rule works.** Moves are proposed at random and accepted by a Metropolis test on a mismatch count
(neighbour pairs holding different rotations) at temperature 0.6, in checkerboard order. The starting knot is built
by hand. The rule is local and stochastic, not a deterministic rewriting rule, and has not been tested for causal
invariance. Record: `src/exgc12_soliton.py`, `src/exgc14_knotgc.py`.

*(Interpretation.)* This is a classical, topological signature. Whether it becomes a quantum sign (a −1 phase under
exchange) is untested and would need a quantum version of the rule.

---

## 6. Motion and inertia: where the programme stalls

### 6.1 Absorption gives drag, not inertia (EX-GC17)
A particle moving through the token gas meets more tokens from ahead. Per unit of inertia the resulting drag is
γ ≈ m₀/6 (within 1–2% at m₀ = 50). Because drag and the gravitational push both come from absorption, the terminal fall
speed is universal (within ±2% of 2E, even where the push alone differs by 7%). The response to acceleration has the
wrong sign for inertia (κ/γ ≈ 0.01–0.07, anti-inertial). So the medium gives Aristotelian, not Newtonian, motion.
*(Interpretation.)* This is the classic objection to Le Sage-type gravity, found here in a concrete model.

### 6.2 Moving the particle as one body
- **EX-GC18 (failed).** Applying the recoil to the knot tears it off the core; pushing W locally squeezes the core. A
  reaction–diffusion pattern has no variable that carries momentum.
- **EX-GC19.** A flow field V, carried by the pattern, lets the particle accelerate as one body; K = 0.2 and 0.4
  particles accelerate within 4%. Cohesion is needed or the core is crushed. The knot cannot keep up.
- **EX-GC20/21.** G as a wave (telegraph equation) and a knot carried with the flow (stays within one lattice step).
- **EX-GC22 (ChatGPT).** Most of the coasting loss in EX-GC19/20 was a defect of the flow model: the transport term
  (V·∇)V was missing. Our earlier "accretion drag" interpretation was withdrawn and errata were added.
- **EX-GC23 (ChatGPT).** With (V·∇)V the particle coasts: the strong particle keeps 97% of its speed over 20 units
  (74% without the term) and its core travels 0.963 of the ideal distance; a weak particle keeps 99.8%. A 2.5% slip
  remains. The inertia is still assigned (τ = 20), and the head-on absorption drag of Section 6.1 was switched off.

### 6.3 The medium's own inertia is far too small (EX-GC24A, ChatGPT)
For a source moving through the wave medium, the reactive (added) mass is positive when damping is weak and negative
when it is strong (consistent with the anti-inertia of EX-GC17). In the weak-damping limit it equals
S²/(6π^{3/2}c²σ), 0.26 at the reference settings. Against the head-on drag b = S, the ratio m/b, which is the coasting
time, is 0.043 for the reference particle and 0.002 for a weak one, against the hundreds of time units needed. The
added mass also scales as S² while gravitational sourcing scales as S, so it could not be a universal mass. Records:
`results/EX-GC24/`, `review/EX-GC24A_claude_review.md`.

---

## 7. The rest-frame obstruction and a frame-free alternative (EX-GC25)

### 7.1 Why consumption-based gravity needs a rest frame
A conserved carrier has a density and a flux that form one four-vector. A uniform background density is that vector
pointing purely in time, which is true in one frame only. A particle moving through it at speed u meets a flux G∞u, and
absorbing that flux costs momentum. This is the drag of Section 6.1, and no choice of medium dynamics removes it.
*(Interpretation, but general.)* Gravity from net absorption of a conserved carrier and Newtonian inertia without a
preferred frame are incompatible.

### 7.2 What survives without a rest frame
- A constant scalar field value is the same in every frame. If G is a scalar field with vacuum value G∞ and a lossless
  wave equation, a particle disturbs it, and the static disturbance obeys Laplace's equation: the 1/r field and the
  1/d² force survive. Like sources of a scalar field attract.
- A particle in uniform motion carries a Lorentz-contracted copy of its field, so there is nothing to drag on.
- If G couples to the particle's energy (the trace of its energy–momentum tensor; Nordström's scalar theory), von
  Laue's theorem (the pressure terms of any stable, localized, static object integrate to zero) makes the gravitational
  charge equal the rest energy, which is the inertial mass. This holds for every particle type and for bound
  composites including their binding energy, which would answer the failures of Section 4.1.

### 7.3 The 1+1 D check
Relativistic φ⁴ kinks stand in for particles. Four kinks with rest energies 8× apart:
- coupled to energy, fell at 0.9987–0.9995 of the predicted universal rate;
- coupled to a charge not tied to energy, spread by a factor of 16, matching Q/M within 0.8%;
- launched at half the wave speed with no potential, lost 2×10⁻⁵ to 2×10⁻⁴ of their speed over 200 time units,
  shrinking with lattice spacing. A kink only 1.8 lattice spacings wide lost 14%: the lattice itself is a rest frame.

Record: `results/EX-GC25/`, figure `paper/figures/exgc25_framefree.png`. The check covers a test body in a fixed
potential only: no self-sourcing, no 3D, no pair.

### 7.4 What the change costs
- **G is no longer conserved or consumed.** The 1/r field becomes a static disturbance, not a flow. This is a new
  model, not a repair of Model Gc.
- **Particles must be lossless.** Friction or a heat bath defines a rest frame, so neither the dissipative
  reaction–diffusion core nor the present Metropolis knot rule qualifies. The knot is the natural particle, but it
  needs a reversible, momentum-conserving rule (and in 3D a Skyrme-type term for stability).
- **Scalar gravity is incomplete.** Nordström's theory gives Newton's law and the equivalence principle but predicts no
  light bending and a perihelion shift of the wrong sign. A full theory needs a tensor part.
- **Discreteness is a rest frame.** Drag returns for particles near the lattice size. Causal invariance is the property
  expected to restore frame independence at large scales in Wolfram-style models; it is untested here.

---

## 8. What is and is not established

**Supported (registered or reproduced across independent checks):**
- 1/r field and 1/d² attraction from a conserved, core-consumed carrier (EX008, EX009, MG195 P1/P2, EX-GC09 on a lattice).

**Shown in exploratory work:**
- Universal free fall in the weak regime under a recoil rule, with a self-energy violation at strong consumption (EX-GC16).
- A stable lattice knot with exactly conserved charge, 720° behaviour and the classical exchange signature (EX-GC12, 13).
- A knot carrying a Gc particle that sources gravity and falls at the right rate (EX-GC14, 15).
- Coasting of a particle given an assigned inertia, once the transport term is included (EX-GC23).
- Absorption drag dominating any medium-generated inertia (EX-GC17, 24A).
- Universal fall and drag-free coasting for energy-coupled scalar gravity, in 1+1 D (EX-GC25).

**Not shown:**
- Newtonian inertia that arises from the model itself.
- Orbits.
- A quantum (rather than classical) exchange sign.
- Causal invariance of any rule used here.
- A frame-free model in 3D.
- Lattice convergence of the particle mass (M differs by about 10% between h = 0.5 and 0.4).

---

## 9. Proposed direction (for the PI's decision)

1. **EX-GC26:** a 3D relativistic lattice knot with a reversible rule and a scalar G field sourced by the knot's
   energy. Tests: 1/d² attraction between two knots; equal fall for two knot sizes; coasting with no drag.
2. **Alternative:** keep Model Gc's consumption picture and look for internal inertia, accepting a finite coasting time.
3. **Later:** a tensor (spin-2) extension; causal-invariance tests of the rule set; a quantum version of the knot rule.

---

## Appendix A. Study index

| ID | Topic | Outcome in one line | Who |
|---|---|---|---|
| EX002–EX007 | Model G gravity (parent repo) | Identity forces zero net G; only screened wells | Claude |
| EX008–EX010 | Model Gc design, pair force, charge | 1/r field, 1/d², mirror charges | Claude |
| MG195 | Inverse-square (confirmatory) | P1, P2 supported; C2 pass; P3 unresolved | Claude, PI closed |
| EX-GC01–04 | Equivalence principle in Gc1 | Non-universal across types and composites | Claude |
| EX-GC05–07 | Rotation / clock | Failed | Claude |
| EX-GC08–11 | Belt trick, tokens, tether, exchange | 720° works; tether gives no exchange sign | Claude |
| EX-GC12–15 | Knot particle | Stable, conserved, exchange sign, carries a falling Gc particle | Claude |
| EX-GC16 | Recoil rule | Universal fall (weak), 1/d² | Claude |
| EX-GC17 | Drag | Drag, not inertia; universal terminal speed | Claude |
| EX-GC18–21 | Recoil on knot; flow field; wave G | Body motion works; coasting loss later traced to a defect | Claude |
| EX-GC22–23 | Flow-model defect; transport term | Particle coasts with assigned inertia | ChatGPT (Claude review) |
| EX-GC24A | Medium-generated inertia | 20–400× too small; scales as S² | ChatGPT (Claude review) |
| EX-GC25 | Frame-free medium | Energy-coupled scalar gravity works in 1+1 D | Claude |

Each study's note is in `results/<ID>/`; plans in `explore/LOG.md` and `spec/`; figures in `paper/figures/`.

## Appendix B. Protocol status
Exploratory studies used the PI's go-ahead in words, an `explore/LOG.md` entry before running and a 2 core-hour cloud
cap. Known deviations are recorded in the log (for example, EX-GC24A ran without a log row and its plan was committed
with its results). Only MG195 is confirmatory. Same-author recomputations are not independent replications.

## References
- Derrick, G. H. (1964). Comments on nonlinear wave equations as models for elementary particles. J. Math. Phys. 5, 1252.
- Finkelstein, D., Rubinstein, J. (1968). Connection between spin, statistics, and kinks. J. Math. Phys. 9, 1762.
- Haisch, B., Rueda, A., Puthoff, H. E. (1994). Inertia as a zero-point-field Lorentz force. Phys. Rev. A 49, 678.
- LaViolette, P. A. Subquantum Kinetics: A Systems Approach to Physics and Cosmology. Starlane Publications (edition to be confirmed).
- von Laue, M. (1911). Zur Dynamik der Relativitätstheorie. Ann. Phys. 35, 524.
- Nordström, G. (1913). Zur Theorie der Gravitation vom Standpunkt des Relativitätsprinzips. Ann. Phys. 42, 533.
- Rueda, A., Haisch, B. (1998). Contribution to inertial mass by reaction of the vacuum to accelerated motion. Found. Phys. 28, 1057 (arXiv:physics/9802031).
- Skyrme, T. H. R. (1961). A non-linear field theory. Proc. R. Soc. A 260, 127.
- Wolfram, S. (2020). A Class of Models with the Potential to Represent Fundamental Physics. Complex Systems 29, 107.
