# Handoff — current state

## EX-GC56 completed (Claude, exploratory, PI CLOSED 2026-11-28)

[Plan](spec/EX-GC56_plan.md) · [Note](results/EX-GC56/EX-GC56_note.md). The shift as a substrate drift. About 2.2 core-h (over the 2 cap, recorded).
- −β·P is exactly a real directed amplitude transfer; it carries every Dirac packet by −β⟨cos k⟩ (7e-5), the same as light. The energy current is mass-dependent (Wilson 1.32–2.22×), so it is excluded by the equivalence principle.
- Momentum-rule pair: κ/2 1.0031, V = 0.15 0.9998, σ = 4 1.0039.
- Light drift P2 narrow fail (0.9942 at k₀ = 0.25), unexplained.
- Code: `src/exgc56.py`, `exgc56_summary.py`, `exgc56_recompute.py`.
- Suggested next (not approved): directed-coincidence amplitude-transfer rule (continuum −β·P, source P, strength 4); ring + matter probe with the synchronised clock.

## EX-GC55 completed (Claude, exploratory, PI CLOSED 2026-11-28)

[Plan](spec/EX-GC55_plan.md) · [Note](results/EX-GC55/EX-GC55_note.md). One interaction for the shift (Codex 54 review). About 0.75 core-h.
- Energy bookkeeping: P/P and reciprocal J_E/J_E conserve; 54's mixed rule does not.
- Pair (n = 49, T = 8): P/P R = 0.9986; reciprocal energy current R = 1.047 (fails). The momentum P is the consistent operator (ADM −β·p).
- Clock: registered gauge test fails (the covariant term is not the cause). Post hoc, preparing the pulse in the clock's locally synchronised frame (`src/exgc55_posthoc2.py`, sign −1) gives gauge residual ≤ 9e-5 and ring wave/readout 0.999–1.005.
- Code: `src/exgc55.py` (covariant loop, je_terms/coupling, evolve with energy), `exgc55_summary.py`, `exgc55_recompute.py`.
- Suggested next (not approved): derive −β·P from directed hop weights; well-depth/speed generalisation with P; ring + matter probe with the synchronised clock.

## EX-GC54 completed (Claude, exploratory, PI CLOSED 2026-11-28)

[Plan](spec/EX-GC54_plan.md) · [Note](results/EX-GC54/EX-GC54_note.md). Source rule for the 53 link shift from Doppler-squared two-site counting (strength 2 × 2 = 4). About 2.7 core-h (PI raised cap to 3).
- Energy-current rule (4κ·½{H, v}): R_dyn(0.2) = 1.019 (registered 1.002 ± 0.008 FAIL; within 2% PASS). Post hoc: in-place strength 3.37 vs 3.85 in uniform N (the current is smaller in the packets' N wells); predicts 1.016.
- Charge-current rule = exactly 2 × energy current (Q = 2H − mβ, {β, α} = 0): R = 0.916, out.
- λ = 2 with 53's law: R = 1.0562, on the line through 52 and 53 (linearity confirmed).
- No-go (N history only): curl 1.2e-4 of needed. Smooth gauge: registered FAIL, post hoc a clock-site β term; propagation is gauge-invariant.
- Code: `src/exgc54.py` (operators Hop/Qop/vop, rule_current, evolve with step-held source), `src/exgc54_summary.py`, `src/exgc54_recompute.py`, `src/exgc54_posthoc.py`.
- Suggested next (not approved): per-local-tick energy current (one extra N factor) in 54B; then the directed-weight coincidence rule.

## EX-GC53 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC53_plan.md) · [Note](results/EX-GC53/EX-GC53_note.md). A signed link shift β (GR's g₀ᵢ), coupled to matter by H = −Σβ·Im(ψ†ψ′) and to light by D_t = ∂_t − β·∇, sourced by the matter current with GR's strength (λ = 4). About 1.6 core-h.
- **Preferred frame removed:** co-moving pair R_dyn = 1.000/0.997/0.989 at V = 0.1/0.2/0.3 (52: 1.115 at 0.2); λ* = 3.91 from 52+53 (GR 4), so the 51 scalar sector is GR-compatible for moving sources.
- **Rotating ring, one-clock light loop:** delay flips with J, zero without rotation, tilt-independent to 0.03%; magnitude → GR dipole as loops grow (registered far-field tests fail: ring radius 10 too large; exact Biot–Savart of the lattice current agrees 1–5%).
- **Still put in by hand:** the shift law and λ = 4 (derivation from firing order did not close).
- **Suggested next (not approved):** a non-symmetric coincidence rule as the shift's mechanism (does it give λ = 4? test with 53E); cleaner far-field ring; later torsion from the 53B spin current, directional strain for tensor waves.

## EX-GC52 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC52_plan.md) · [Note](results/EX-GC52/EX-GC52_note.md). Moving source in the 51 model (PI go-ahead 19:54 UTC). About 1.0 core-h.
- **Field of a moving mass is right:** Lorentz-contracted with charge Q(V), g₀₀ = boosted static (amplitude 1.000 after ±y symmetrisation).
- **Preferred frame confirmed:** co-moving pair R(V) = γṗ(V)/ṗ(0) = 1.02/1.11/1.28 (V = 0.1/0.2/0.3; symmetrised 1.03/1.13/1.33 = lattice prediction to 0.4%); full dynamic R_dyn(0.2) = 1.115. Lorentz invariance needs 1. Missing g₀ᵢ (PPN α₁ ≈ −8, α₂ ≈ −1).
- **Side finding:** moving spin-up Dirac packet's charge centre shifts −0.18 sites in y at V = 0.2 (reverses with spin; Møller-shift size). This made the P3 1% match and P5 fail.
- **Suggested next (not approved):** a minimal direction-carrying addition (torsion / handed coincidences sourced by matter current, per PI's idea) and rerun 52D (0.4 core-h) to see if R → 1.

## EX-GC51 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC51_plan.md) · [Note](results/EX-GC51/EX-GC51_note.md). Codex 48/49/50 review recommendation 3: one declared model (Wilson–Dirac matter with coincidence weights, N carried as light's slowness, sourced by q = ∂H/∂N, light by 37B), run end to end. About 6.5 core-h (PI raised cap twice).
- **Passes:** G_eff = 1 to 3e-5; field and light arrive together to 0.1%; reciprocity |F₁+F₂|/|F₁| = 0.014; ray Shapiro ratio 2.02–2.04.
- **Registered fails:** c_N = 1.08 (light reads the same; likely the Dirichlet-box target); β separation as coded 1.06 (used κ not κQ; post hoc 1.000); wave Shapiro 2.076 at b = 8 (second-order U); all P3 clauses (comparator omitted ⟨E_k E_kN⟩∂N; σ = 3/m probes differ in size, so r = 0 equal fall was untestable).
- **Post hoc:** with the velocity term, fall matches the lattice dispersion within 2–6%; point limit (1 + rm)×GR.
- **Suggested next (not approved):** moving source (does a scalar give g₀ᵢ effects?); radiation from an orbiting pair; clean universality rerun with equal σ; an action for the coupled rule.

## EX-GC50 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC50_plan.md) · [Note](results/EX-GC50/EX-GC50_note.md). Idea 4 from the GC40 recap: a natural reason for 42's vacuum subtraction.
- **Clock sector: yes, exact.** Under 44's homogeneous weights the vacuum has no lapse mass (Euler) and χ_NN ≤ 0 (concavity); network q² stiffness −0.0025/−0.010 (vs 43's +1.00); Dirac sea only q⁴.
- **Spatial sector: no.** χ_Na(0) = −ε_link and χ_aa(0) (equal to 43's 37B χ(0) exactly) remain: the counterterm of 42–47 lives there. Lattice cosmological-constant problem.
- **Only exact boson–fermion partners cancel it** (Dirac sea + partner, zero to 1e-16); unbroken lattice SUSY. 47's G_eff would be ~1% larger; conclusions unchanged.
- Fails recorded: Dirac χ₀ clause at L = 32 (q⁴ behaviour; post hoc L = 96 fine).
- **Suggested next (not approved):** a spatial-factor coupling homogeneous like the clock's (does Q_a and light's factor 2 survive?), or softly broken partners.

## EX-GC49 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC49_plan.md) · [Note](results/EX-GC49/EX-GC49_note.md). Idea 1 from the GC40 recap: firing-rate changes carried by the same two-site coincidences as light.
- **The rule derives a wave equation for N and shares the light cone** (light and δN Shapiro delays agree to 0.05% at k₀ = 1).
- **But passing on N itself gives PPN β = −1/2** (N³ harmonic): perihelion 1.5× GR, ruled out. Of four transported quantities only the light slowness 1/N² gives β = 1 (exponential metric).
- Scored lattice fit missed ±0.05 for 3 of 4 from a linear-in-strength lattice anisotropy artefact; post hoc separation recovers the continuum β to 0.006.
- **Suggested next (not approved):** a reason in the event picture for passing on delay (slowness) rather than rate; large-amplitude behaviour (no horizon in the exponential metric).

## EX-GC48 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC48_plan.md) · [Note](results/EX-GC48/EX-GC48_note.md). Idea 7 from the GC40 recap (Codex GC47R test): a Wilson–Dirac lattice fermion under the 44 rules.
- **One fermion is relativistic matter:** 2 t_pair E_rest/c² = 1 + r m (exact), with c measured from massless packets; dust at rest (Q_a/E ∝ σ⁻²); moving stress ∝ v². The whole Wilson term must stretch with a (W1); otherwise Q_a/E = −3r/m.
- **47 calibration:** with the network's own fastest signal (0.174), the 47 pairs have 2 t_pair E/v² ≈ 0.004. The inertia failure stands.
- **A contact-bound composite is not:** inertia 1.27–3.6 times its energy/c² (corrected 2026-11-27; R_c = 0.28–0.79 is the reciprocal), binding stress η = 40–318. Contact binding is cutoff-scale (T, |V| 2–3 ≫ m); instantaneous potential energy carries no inertia or stress.
- **Suggested next (not approved):** binding by exchange of a propagating lattice field, to test whether field momentum/stress restore R_c = 1, η = 0.
- Fails recorded: P2 at m = 0.2 (asymptotic formula scored too early), P3b one point, P3d at m = 0.2, C1 Q_a finite-difference noise.

## EX-GC47R completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC47R_plan.md) · [Note](results/EX-GC47R/EX-GC47R_note.md). Corrections to 47 after the Codex review. The release packet is now periodic about its own centre, and the form factors are Hellmann–Feynman consistent (projector-derivative term). 47's conclusions stand: acceleration = −2 t_pair ∂V/∂r within 4.4%, a/a_N ≈ 1.2 × 10⁻⁴, force = point-charge force to 4e-4 and → 1/r² with box size. The energy-drift gate failed (split-step, 2.3e-7). Next per Codex (not approved): a massive band tested against 2 t_pair E_rest = c_signal², with c_signal measured independently.

## EX-GC47 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC47_plan.md) · [Note](results/EX-GC47/EX-GC47_note.md). Implements the Codex GC46 review's three-stage two-body test.
- **Two bound heavy pairs attract** through the 44 fields exactly as the point-charge formula −G(E₁E₂ + E₁Q_a2 + Q_a1E₂) predicts (to 0.14% at L = 32), reciprocally, and the attraction moves them (acceleration = −2 t_pair ∂V/∂r to ≤ 2.6%, validated against the exact two-body problem at L = 6 to 1.9%).
- **But they fall at about 1.2 × 10⁻⁴ of the Newtonian rate.** A pair's inertia 1/(2 t_pair) ≈ 8,000 is not tied to its energy ≈ 1. The on-site rest energy that made matter dust-like (45/46) is what breaks m_i = E.
- **Two clauses fail:** F_a(0) = Q_a (0.8%: projected model breaks strain Hellmann–Feynman), and B1 vs B2 free fall differs by 5–9%, not > 20%.
- **Suggested next (not approved):** a Dirac-type (Wilson) lattice fermion with relativistic dispersion. It would give m_i = E/c² and Q_a = 0 at rest without a small parameter; test its bound pairs for m_i = E and zero binding stress.

## EX-GC46 completed (Claude, exploratory, PI CLOSED 2026-11-27)

[Plan](spec/EX-GC46_plan.md) · [Note](results/EX-GC46/EX-GC46_note.md). Follows the Codex review of 45 (tests 1 and 2). 43, 44 and 45 were PI CLOSED 2026-11-26 (Richard's go-ahead).
- **Four upper bands per carrier** (rank-1 contact in a larger basis) keep the binding-stress result: η = δQ_a/E_b is −2.3 to −1.1 at strong hybridisation and +1.1 to +4.9 at weak. At the middle point it is −0.033 (Codex: −0.039). The projection error in E_b is 13–18% at strong hybridisation and 0.3% at weak.
- **Small-g limit:** the pair's total spatial charge is 0.2–0.3% of its energy at g = 0.0125, and its mobility falls as g². η tends to a finite positive value (0.4–2.7). The sign flip belongs to g, the size to m₀.
- **Suggested next (not approved):** Codex test 3, a finite-range interaction with a pre-fixed strain law, to see whether η can be made to vanish; or a lattice dilation that scales interaction and hopping together.

## EX-GC45 completed (Claude, exploratory, PI CLOSED 2026-11-26)

[Plan](spec/EX-GC45_plan.md) · [Note](results/EX-GC45/EX-GC45_note.md). Follows the Codex review of 44 (items 1 and 2; item 3, a two-field optics test, deferred).
- **Heavy on-site matter is nearly dust.** A rest mode on a site, reached by hybridisation g, has Q_a/ΔE ≈ 0.008 at g = 0.025, falling as g². But Q_a stays 0.3–0.7 of its band width: dust only in the heavy (non-relativistic) limit.
- **Binding energy does not gravitate as energy.** A bound pair has a lattice Nordtvedt ratio η = δQ_a/E_b of −2.7 to +4.9 (GR: 0), with a sign set by the hybridisation. No interaction strain exponent cancels it universally (−0.6 to +1.4).
- **C5 (projection error) failed at its registered point**, where η is near its sign change; elsewhere the error is small enough that the sign and size of η stand.
- **44's exciton at its true minimum** stays radiation-like (Q_a/ΔE 1.40–2.17).
- **Suggested next (not approved):** an interaction that acts across links with the same strain factor a as hopping (a lattice dilation), to test whether a von Laue identity then holds; or Codex's two-field optics test.

## EX-GC44 completed (Claude, exploratory, PI CLOSED 2026-11-26)

[Plan](spec/EX-GC44_plan.md) · [Note](results/EX-GC44/EX-GC44_note.md). Follows the Codex review of 43: separate the clock (lapse N, on everything) from space (link factor a).
- **Lapse-homogeneous weights:** mass equals energy for every state, and the redshift test passes. But light bends at only the Newtonian value.
- **A spatial link factor:** gives a second charge, the strain stress Q_a = ΔE_link (the lattice 3p). 43's charge is Q_N + Q_a.
- **This network's matter is radiation-like:** Q_a ≈ ΔE even at rest.
- **A mobile exciton binds,** but its stress cancels only if the attraction grows as links stretch (w* < 0). No ordinary interaction does that.
- **Suggested next (not approved):** matter whose rest energy is on-site (Q_a = 0 at rest). Test whether its composites are stress-free.

## EX-GC43 completed (Claude, exploratory, PI CLOSED 2026-11-26)

[Plan](spec/EX-GC43_plan.md) · [Note](results/EX-GC43/EX-GC43_note.md). This follows the Codex review of 42: the lapse couples to the conjugate source ∂H/∂N, derived from one functional.
- **Mass is not energy.** The charge is Q = ΔE + ΔE_link. Free packets attract at 1.8–2.25× their energy (a pair up to 3.2×), and hubs at 0.10–0.41×. No state has Q = ΔE, so the equivalence principle fails at order 1. It is the same defect as 42A's redshift failure.
- **The vacuum term needs a counterterm.** χ(0) ≈ −6 makes the raw constraint Helmholtz and boundary-set. The local counterterm gives 1/r, with G screened by c₂χ ≈ 1.
- **Inference.** GR gets Q = E for bound bodies from von Laue (stresses integrate to zero). This lattice has no dilations and no self-bound states, so nothing forces it.
- **Suggested next (not approved):**
  - an interacting network with bound composites (does the excess link energy vanish for bound states?);
  - or lapse weights designed so that the conjugate charge equals the energy for stationary states.
- **Recommendation:** stop scaling up the Gaussian pair network.

**2026-11-26 closures.** Richard: "go ahead and close anything that is open" (11:05 America/Denver). Every Model Gc exploration still open is now PI CLOSED: EX-GC01–26B3 (rows with no closure recorded), 35Q, 35U–X, 37A-S2, 37C, 37D, 38A, 38B, 39A–C, 40A–C, 41A/B and 42. No result was changed. Nothing is open now; the next study needs a new plan and go-ahead.

## EX-GC42 completed (Claude, exploratory, PI CLOSED 2026-11-26)

PI go-ahead 2026-11-26 17:27 UTC. [Plan](spec/EX-GC42_plan.md) · [Note](results/EX-GC42/EX-GC42_note.md). The 40A pair network is coupled to a lapse N by 37B weights (links N_vN_w, sites N_v), and N is fixed by lattice Poisson sourced by energy density.
- **With the local vacuum energy subtracted (variant A), a hub sources a universal 1/r.** It holds for 1 or 7 sites and every Δλ, λ₀ and G: the coefficient is the same to 3e-5, the mass equals ΔE₀, and G_eff = G/(1 + 4πG c₂) with c₂ ≈ 0.49 (mostly the half-link energy bookkeeping).
- **If all energy gravitates (variant B), there is no Newton law.** dρ_vac/dN ≈ −6 makes the constraint a Helmholtz equation, which is box-resonant with no unique decaying solution. This is the cosmological-constant problem in miniature.
- **A hub's mass is negative** (sharing lowers the energy), so it raises N nearby.
- **The EP fails under 37B weights** (ΔE₀/N changes 11–17% over N = 0.9–1.1).
- Suggested next questions (not approved): a natural reason for the vacuum subtraction; whether a positive-energy excitation sources the same G_eff; an N-independent hub.

## EX-GC41 and EX-GC40C completed (Claude, exploratory, PI CLOSED 2026-11-26)

These follow the Codex review of 40A/B. Plans: [EX-GC41](spec/EX-GC41_plan.md), [EX-GC40C](spec/EX-GC40C_plan.md).
- **[41](results/EX-GC41/EX-GC41_note.md): delay works, bending does not.**
  - On a random network, the front delay matches the continuum lapse calculation to 2% (the Shapiro delay, recovered absolutely).
  - Bending from the front slope is swamped by first-passage roughness.
  - With shared site clocks, the median-of-histories two-site delay is 16% above the mean-field value. So GR's factor 2 is a mean-field property of the rule.
- **[40C](results/EX-GC40C/EX-GC40C_note.md): a hub is a local defect.**
  - It binds states and two hubs attract, but all responses decay exponentially (about ξ_C/2), with no 1/r.
  - So free pair networks cannot make long-range gravity.
  - Suggested next step (not approved): couple a pair network to a lapse field sourced by its energy density, and ask whether a hub sources a universal 1/r lapse.

## EX-GC40A and 40B completed (Claude, exploratory, PI CLOSED 2026-11-26); 40C outlined, not approved

Richard shared the Chain Theory compendium (2026-11-26). Its gravity is unfinished: its own field equation fails at all 60 cells it tests, and inertia is not treated. We took only its axiom, a vacuum made of maximally entangled pairs. Plan: [spec/EX-GC40_plan.md](spec/EX-GC40_plan.md).
- [40A](results/EX-GC40A/EX-GC40A_note.md):
  - Pure pairs do not respond to matter at first order.
  - Once pairs share entanglement with their neighbours (λ > 0), area as "cut links" works: rewiring changes S by s_link per link, within 3%.
  - The matter response is attractive but grows like ℓ^2.4–2.9, so G_eq is not constant (it falls 2.3–3× over ℓ = 3–6).
  - Entanglement equilibrium has now failed to give a constant G on three substrates.
- [40B](results/EX-GC40B/EX-GC40B_note.md):
  - 37B's counting rule works on a random network with no lengths: v ∝ N^1.996, and the Shapiro ratio is 2.004.
  - The bending readout was too noisy, and its target was not 2 for fixed endpoints.
  - The isotropy gate was below the noise floor.
- 40C (localised shared-entanglement "hub" as a mass) is outlined in the plan. It needs its own predictions and a go-ahead.

## EX-GC39A, 39B and 39C completed (Claude, exploratory, PI CLOSED 2026-11-26)

Richard's question (2026-11-26): does matter set the lapse, or does the lapse set how matter emerges? Design: [spec/EX-GC39_design.md](spec/EX-GC39_design.md).
- [39B](results/EX-GC39B/EX-GC39B_note.md): the lapse can shape matter invisibly (same proper size and clock; geodesic fall at 0.997) only if matter is wave-built from the same events as light. A Skyrmion needs its quartic term at weight N⁶. Model Gc's reaction–diffusion particles fail (proper size ∝ N^-1/2). A self-sourced lapse well holds its own field (the Schrödinger–Newton family), up to GμM ≈ 0.55, with no mass selected.
- [39C](results/EX-GC39C/EX-GC39C_note.md): a free wave lapse radiates 1/6 of GR's quadrupole power, so binary pulsars exclude a free N. The lapse must be constrained by matter, as in GR.
- [39A](results/EX-GC39A/EX-GC39A_note.md): entanglement equilibrium on a fixed lattice gives a G, but its size and even its sign depend on how "fixed volume" is defined. The responses are boundary effects. For geometry to come from entanglement, the counts of sites and links must themselves be the geometry.
- Next idea (not approved): a graph or causal-set substrate where area = number of cut links, with the 39A test repeated there.

## EX-GC38A and 38B completed (Claude, exploratory, PI CLOSED 2026-11-26)

Richard's question (2026-11-26): is what makes matter gravitate the mechanism of entanglement?
- [38A](results/EX-GC38A/EX-GC38A_note.md): vacuum entanglement per unit area vs 37A's induced 1/(4G), same lattice scalar. The universal m² ln m terms agree exactly (1/24π, for both cut orientations and both graviton polarisations). The cutoff-dominated parts do not: the ratio is 2.0–4.8. The plan's continuum value had a factor-2 slip, recorded as a failed clause.
- [38B](results/EX-GC38B/EX-GC38B_note.md): the half-chain entanglement Hamiltonian is the boost generator, to 0.02–0.8% (lattice Bisognano–Wichmann). With a lapse bump, the local entanglement temperature moves the Tolman way but only about 60% as much as the optical-length prediction. That is a P2 FAIL, and it is not caused by the mass. Open.
- Next idea (not approved): a Jacobson design note, i.e. use δS = δ⟨K⟩ on small cuts to derive the equation for the lapse N.

## EX-GC37C, 37D and 37A-S2 completed (Claude, exploratory, PI CLOSED 2026-11-26)

- [37C](results/EX-GC37C/EX-GC37C_note.md): moving source. A diffusing, consumed carrier gives a wake and self-drag, i.e. a preferred frame. A wave field with speed 1, sourced by the same consumption, is fore/aft symmetric, Lorentz-contracted (γ within 1.3% to v = 0.4) and drag-free. The 37B mechanism needs N to be a relativistic wave field.
- [37D](results/EX-GC37D/EX-GC37D_note.md): realising coincidences randomly, event by event, on a wave pumps energy and destroys coherence. A coherent wave has to be a coarse-grained mode of many coincidences.
- [37A S2](results/EX-GC37A/EX-GC37A_note.md): static scalar sector of induced gravity on Y1. γ = 0.236 (gauge h_L = 0), large Ward violations, but the lapse stays non-dynamical. Fails worse than the spin-2 sector.
- Next idea (not approved): a rule where firing-rate changes propagate by the same two-site coincidences as light, so N obeys a wave equation at the light speed.

## EX-GC37A and EX-GC37B completed (Claude, exploratory, PI CLOSED 2026-11-26)

Both test a *mechanism* for gravity rather than putting GR's equation in by hand (Richard's restated goal, 2026-11-25).
- [37B](results/EX-GC37B/EX-GC37B_note.md): sites fire at a local rate N. Clocks are one-site events (rate N). Propagation needs two sites to fire together (speed N²). That gives the weak-field metric with γ = 1: token delay ratio 2.00–2.08, wave bending ratio 2.04, no spatial metric put in. A consumed, diffusing carrier sets N = 1 − GM/r. Open: why propagation is a coincidence, second order, and the moving-source / preferred-frame test (B3).
- [37A](results/EX-GC37A/EX-GC37A_note.md): Sakharov induced gravity on Y1, spin-2 sector. The right sign, but gravitons at 0.78–0.85 c, the two polarisations 2.1× apart (corrected after ChatGPT review), and a huge mass term. The Lorentz-invariant continuum control gives c_g = 1. So an induced-stiffness substrate needs a locally Lorentz-symmetric cutoff. S2 (γ, modes) not run.

## EX-GC36C completed (Claude, exploratory, PI CLOSED 2026-11-25)

[Results](results/EX-GC36C/EX-GC36C_note.md). Rule 66427's order-free particle does not persist in its own debris. The median lifetime is 3 events, and only 5–16% of schedules are alive at 500, every one of them as a proliferating swarm of particle-shaped supports. The 55744 control lives longer and shows the same pattern. Combined with 36A/36B, the CRR 80-rule catalogue has no order-free persistent particle and no particle whose forks rejoin. The 36B vacuum join count is a lower bound for version L (ChatGPT review).

## EX-GC36B completed (Claude, exploratory, PI CLOSED 2026-11-25)

Follows [ChatGPT's 36A review](review/EX-GC36A_chatgpt_review.md) ([results](results/EX-GC36B/EX-GC36B_note.md)). Only CRR rule 66427 of 80 keeps a closed particle without the vertex-order tie-break. The odd 36A forks come from duplicate-output suppression. A state-frontier test, calibrated on known delayed joins, confirms the published Wolfram rule but finds no rejoining for the 55744 or 66427 particle forks, nor for 27 of 28 vacuum forks, within 5–6 events per side.

## EX-GC36A completed (Claude, exploratory, PI CLOSED 2026-11-25)

Causal-invariance audit of CRR Rule 55744 ([design](docs/EX-GC36_rule55744_causal_invariance_design.md), [results](results/EX-GC36A/EX-GC36A_note.md)). The calibration gate was not met: the published overlapping rule joins only 1 of 3 forks on an asymmetric seed within 3 steps. No overlapping 55744 fork joined, in the vacuum (0/28) or at the particle's B fork. CRR's particle closure depends on the u < v birth-order tie-break. Version C of Stage 2a and Stage 2b were not completed.

## GC35U–X completed

Ran items 1–4 of the [GC35R–T review](review/EX-GC35RST_claude_review.md) in the requested order. These studies are exploratory and pending PI closure.

- **GC35U:** restricting the initial anchors to producible symbols recovered the two ternary string passes. Adding left or right transport broke their plain causal graphs. The string line is complete. [Results](results/EX-GC35U/EX-GC35U_note.md)
- **GC35V:** 236 of 297 typed rules had shared-edge forks, and none found an anchored join within the registered bound. Sixteen first witnesses had terminal arms. No typed rule reached the conditional extension. [Results](results/EX-GC35V/EX-GC35V_note.md)
- **GC35W:** changing P's position can make an eligible absorption disappear. That produces different aggregate ledgers and causal graphs with both recoil thresholds. [Results](results/EX-GC35W/EX-GC35W_note.md)
- **GC35X:** gating the probe on a local background event enforced one site per gated tick by construction. Doubling fixed-site rate changed the probe's speed by only about 4%. [Results](results/EX-GC35X/EX-GC35X_note.md)

The inflow, attraction and continuum knot lines remain paused.

## GC35R–T completed

- **GC35R** screened 734 nonexpanding singleton or two-rule systems. Four passed the output-border and plain graph screen through length seven. Zero passed with both input and output ports anchored. [Results](results/EX-GC35R/EX-GC35R_note.md)
- **GC35S** supersedes the earlier exclusion clock stop under the amended scope. Blocking is an interaction without shared consumed sites. All common depth cuts agreed across schedules. Forward density calibration fit within 1.5%, while the step showed forward/reverse rate differences and the independent probe clocks differed by 42% in inferred speed. [Results](results/EX-GC35S/EX-GC35S_note.md)
- **GC35T** preserved parallel absorption by keeping P as a persistent vertex and consuming independent carrier relations. All 124 configurations passed. A shared consumed mass counter failed the anchored and routed graph comparisons. [Results](results/EX-GC35T/EX-GC35T_note.md)

No configuration combines qualified anchored competing histories with an absorber transport rule ready for the gated clock, inflow or attraction tests. GC35R–T are exploratory. **PI CLOSED 2026-11-25.**

## Emergent clocks — Step 2 held at the entry gate

The [Step 2 design](docs/EX-GC35_emergent_clocks_design.md) requires an interacting rule that passes the Step 1 graph screen and directs the run to stop when none qualifies. GC35Q passed only the one-species exclusion control, which has no competing matches. Pass-through and both absorber variants have plain causal-graph counterexamples. No clock run was started and no density–rate or probe-speed result was generated. [GC35Q results](results/EX-GC35Q/EX-GC35Q_note.md)

## GC35Q completed

GC35Q ran Step 1 of the [GC35P review](review/EX-GC35P_claude_review.md). The pass-through overlap R0L has equivalent plain graphs, while the four-site context 0R0L has unequal causal depths despite identical terminal words. Both absorber variants have larger-context counterexamples. Only the non-overlapping exclusion control passed. GC35Q is exploratory and pending PI closure. [Results](results/EX-GC35Q/EX-GC35Q_note.md)

## GC35L–P closure

PI CLOSED 2026-11-24 (Richard: “close 35L through 35P”). The five outcomes are recorded in the [experiment log](explore/LOG.md).

## GC35L–O review follow-up

- **GC35P** screened five carrier rules. Independent directed hops conserve carriers and preserve both tested graphs. In the reduced reciprocal clock loop all 16 randomized schedules matched, including fields, event times and labelled dependencies. Each rule with contested inputs or queue batching had a graph counterexample. Finding an interacting rule that passes is the open carrier problem. [Results](results/EX-GC35P/EX-GC35P_note.md)
- **Knot-on-Y1 continuum work paused.** Energy-only sourcing is compatible to roundoff. The full-stress residual falls roughly as h² over the measured spacings and remains nonzero. The local damping test did not restore compatibility with an unconserved source. The exact identity failed for the tested hard and soft couplings. Revisit this line if a carrier result calls for a continuum check. [Damping](results/EX-GC35M/EX-GC35M_note.md) · [Conservation note](docs/EX-GC35_local_matter_conservation.md)
- **Lens campaign complete.** The full-to-lapse factor remained 2.007–2.009. The dispersion-corrected absolute gap stayed around 4.4% at the last two spacings. The finite-beam versus straight-ray comparator is a plausible explanation. The existing runs do not isolate it. [Finest lens test](results/EX-GC35N/EX-GC35N_note.md)
- **White paper v0.5 frozen through GC35O** in revision 7. GC35L–P are exploratory and PI closed. [White paper](paper/Model_Gc_white_paper_v0.5.md)

## GC35A–K closure and GC35L–O follow-ups

PI CLOSED 2026-11-24 (Richard: “close 35A through 35K”). All eleven note statuses and log entries were updated with their original outcomes preserved. The [white paper v0.5 revision 6](paper/Model_Gc_white_paper_v0.5.md) includes GC35J/K.

- **GC35L** found equal ledgers and firing counts in 494 configurations, with inequivalent token and FIFO site-history graphs. [Results](results/EX-GC35L/EX-GC35L_note.md)
- **GC35M** lowered the fine scalar residual to 0.0449 at t=20, above its gate. The reconstructed tidal error was 27.6%. The fixed full-stress knot sources remained incompatible while the energy-only controls stayed at roundoff. [Results](results/EX-GC35M/EX-GC35M_note.md)
- **GC35N** extended the lens to h=.08. The full/lapse bend ratio was 2.007 and the dispersion-corrected full-angle gap remained about 4.4%. [Results](results/EX-GC35N/EX-GC35N_note.md)
- **GC35O** demonstrated two-way packet and token lapse feedback in a one-dimensional static-response reduction of Y1. Field schedules agreed under fixed version semantics, while contested token branches gave different causal graphs despite identical counts and fields. [Results](results/EX-GC35O/EX-GC35O_note.md)

The exact-Ward knot search is paused. The [literature note](docs/EX-GC35_local_matter_conservation.md) records why the known nonlinear Klein–Gordon result has a narrower scope than a general barrier. GC35L–O are exploratory and pending PI closure.

## Follow-up to the GC35G–I review

- **GC35G** is recorded as the working projected scalar benchmark for the tested window. The refinement campaign is complete. The recipe marks the Fourier projection as non-local and the h=.0625 value as an extrapolation. PI closed 2026-11-24. [Recipe](docs/EX-GC35_projected_Y1_recipe.md)
- **GC35J** failed the direct soft-potential Ward gate. The quartic radial potential does not obey the endpoint-averaged link chain rule. The gated relaxation and moving-knot runs were not started. [Results](results/EX-GC35J/EX-GC35J_note.md)
- **GC35K** used a static matter-sourced Y1 metric with local proper-time events. Replacing the linear predictor with causal quadratic version reconstruction recovered the uniform-clock angle control. Fine group bends were 3.135° full and 1.560° lapse-only, versus 2.834° for the matched finite-path full prediction. [Results](results/EX-GC35K/EX-GC35K_note.md)
- The [competing-event design](docs/EX-GC35_competing_event_causal_invariance.md) defines an exact graph-equivalence comparison and a weaker ledger comparison for a future run. [White paper v0.5](paper/Model_Gc_white_paper_v0.5.md) now covers the research record through GC35K.

## GC35G–I completed

Items 1–3 of the [GC35E/F review](review/EX-GC35EF_claude_review.md) were run in the PI's requested sequence. All three are exploratory and PI closed 2026-11-24.

- **GC35G** confirmed near-fourfold refinement of the reconstructed projected tide at h=.125. The late error was 2.25% at t=20, above the 2% gate. Both constraint gates passed. The [recipe](docs/EX-GC35_projected_Y1_recipe.md) and [v0.5 draft revision 4](paper/Model_Gc_white_paper_v0.5.md) were committed before H and I began. [Results](results/EX-GC35G/EX-GC35G_note.md)
- **GC35H** reduced the static h=.2 full stress residuals to 0.962% and 0.679%. The exact Ward gate failed. A stationary rank-one configuration isolates a constraint reaction term in the native sigma source. The moving-knot run remains gated on compatibility. [Results](results/EX-GC35H/EX-GC35H_note.md)
- **GC35I** passed the fixed-version schedule and continuum comparisons in the prescribed lapse step. All 18 randomized schedules agreed with the coordinate-ordered control, with fast/slow event density 1.5. At h=.1 the proper-frequency ratio was 0.667113 and the outgoing angle 32.851°. [Results](results/EX-GC35I/EX-GC35I_note.md)

## GC35E/F completed

- **GC35E** used a tensor-weighted projection after each step of the GC35C reciprocal run. Post-projection residuals were at roundoff through t=20. The fine pre-projection peak of 1.86e-3 and reconstructed tidal error of 8.94% at t=20 failed the accuracy gates. The TT evolution was unchanged by the projection to numerical precision. [Results](results/EX-GC35E/EX-GC35E_note.md)
- **GC35F** reproduced the one-way linear Y1 scalar-metric stencil on a fixed versioned local dependency graph without a timestep barrier. All 48 legal schedules agreed bit for bit and the stale-read negative control was schedule dependent. Common-slice constraint errors fell by 4.01× on halving the timestep. This establishes order independence within the prescribed graph. [Results](results/EX-GC35F/EX-GC35F_note.md)

These were items 1 and 3 of the [GC35C/D review](review/EX-GC35CD_claude_review.md). Both studies are exploratory and PI closed 2026-11-24.

## GC35C/D completed

- **GC35C** extended the compensated scalar benchmark to t=20. The fine residual exceeded 10⁻³ at t=5.75 and reached 0.368 at t=20. Neither the timestep change nor removing the compensator corrections accounted for the growth. The linear Y1 vacuum generator had no matching exponential mode. This limits the current reciprocal scalar construction to its qualified short window. [Results](results/EX-GC35C/EX-GC35C_note.md)
- **GC35D** met the tighter rest-force tolerance in all eight cases. Matched relaxation and box controls hardly changed the residual at h=.2. At h=.15 the full-stress residuals were 0.867% and 0.698%, above the 0.1% gate. [Results](results/EX-GC35D/EX-GC35D_note.md)

The moving-knot study was conditional on both gates and was not started. GC35C and GC35D are exploratory and PI closed 2026-11-24.

## GC35A/B and white paper v0.5 revision 3

The PI accepted the approximate staggered Y1 model as the current reference and paused exact nonlinear correction searches. [Symmetry and formulation note](docs/EX-GC34_symmetry_and_direction.md) · [White paper draft through GC34Z](paper/Model_Gc_white_paper_v0.5.md)

- **GC35A** completed. A0–A2 passed in the finite compensated periodic scalar run. Halving spacing reduced the final next-order constraint defect by 3.31×, with negligible change from timestep halving. Fine-grid tidal errors were 1.35% and 5.27% at the two readouts. The defect grows within the tested window and its long-time behavior is open. [Results](results/EX-GC35A/EX-GC35A_note.md)
- **GC35B** completed. The registered knot comparison failed. At post hoc h=.2 both knots survived and gave controlled energy-only field and initial feedback results. The full energy-plus-stress source missed the compatibility gate for both. [Results](results/EX-GC35B/EX-GC35B_note.md)

Both GC35 studies are exploratory and PI closed 2026-11-24.

## EX-GC34S, X, Y and Z closure

PI CLOSED 2026-11-24 (Richard: “I sign off on all of the open items”). GC34S ruled out direct reciprocal sourcing with the linear kernel. GC34X–Z ruled out the tested cubic interaction families at the registered tolerance. [GC34S](results/EX-GC34S/EX-GC34S_note.md) · [GC34X](results/EX-GC34X/EX-GC34X_note.md) · [GC34Y](results/EX-GC34Y/EX-GC34Y_note.md) · [GC34Z](results/EX-GC34Z/EX-GC34Z_note.md)

## EX-GC34V and EX-GC34W closure

PI CLOSED 2026-11-24 (Richard: “go ahead and close V and W”). GC34V's tested quadratic Y1 source failed finite-spacing symmetry and conservation. GC34W's variational cubic action gave a symmetric source but failed the vacuum conservation gate. Neither candidate went to reciprocal evolution. [GC34V results](results/EX-GC34V/EX-GC34V_note.md) · [GC34W results](results/EX-GC34W/EX-GC34W_note.md)

## EX-GC34U closure

PI CLOSED 2026-11-24 (Richard: “Close GC34U.”). The isolated spherical continuum reciprocal benchmark passed U0–U2 on the corrected implementation. [Results](results/EX-GC34U/EX-GC34U_note.md)

## EX-GC34R closure

PI CLOSED 2026-11-24 (Richard: “Go ahead and close 34R.”). All three registered gates passed for the scalar packet’s leading response to its own weak field. The gravity source was held at zeroth order. [Results](results/EX-GC34R/EX-GC34R_note.md)

## EX-GC34T closure

PI CLOSED 2026-11-24 (Richard: “go ahead and close 34T.”). Both registered checks passed. GC34T verified the continuum next-order conservation identity and the expanding homogeneous benchmark. The discrete Y1 identity and reciprocal packet evolution remain open. [Results](results/EX-GC34T/EX-GC34T_note.md)

## EX-GC34Q closure

PI CLOSED 2026-11-23 (Richard: “we can close GC34Q”). GC34Q passed Q0–Q3. The differential tidal curvature from added scalar kinetic energy matched linearized GR within 0.771% on the finer grid, with second-order convergence. This was a one-way weak-field test with the scalar evolving on the flat grid. [Results](results/EX-GC34Q/EX-GC34Q_note.md).

## Current state (2026-11-22, updated after PI closure of EX-GC31A–33A)
**No study is running.** EX-GC26 (26A–F, 26E3), EX-GC27A/B, EX-GC28A–H, EX-GC29A and EX-GC30A/B/B2/B3 PI CLOSED (2026-11-22; 30-series review review/EX-GC30A-B3_claude_review.md). EX-GC29B deferred; white paper v0.4 covers through EX-GC26.
**EX-GC30:** analytic screen, BEC tuning failure under density, Weyl winding conserved but metrics split except symmetry-restricted controls; weak RG common-speed attraction insufficient at finite scales. PI CLOSED. Nonzero interacting fixed point remains reserve.
**EX-GC31A (ChatGPT):** analytical test complete; review and PI closure pending. Exact asynchronous deterministic port-rule examples prove schedule-invariant event DAG and shared structural reach do not enforce common operational geometry. A: sector fronts 1 vs 1/2; B: identical extremal microscopic fronts 1 with gapless longwave speeds 0.866 vs 0.5 under exchange symmetry. Positive control imposes common principal symbol rather than deriving it. Parent ModelG MG189 causal completion remains unresolved; no parent changes. See results/EX-GC31A/EX-GC31A_note.md and spec/EX-GC31A_plan.md. Negligible compute, paper only.
**EX-GC31B1 (ChatGPT):** 80 archived Stage W rules screened on 179 constructed critical overlaps each using documented reconstructed grammar (original engine ZIP local-only). 80/80 disjoint controls pass; 0/14,320 rooted two-event diamonds; no terminal mismatch; coarse dependency-only matches exist in every case. H1 remains unresolved globally. B2/B3 not opened. results/EX-GC31B1/EX-GC31B1_note.md. Review/PI closure pending.
**EX-GC31B1B (ChatGPT):** immutable relational pair join, exact fair-completion and relabel-covariance proof, 120/120 exact cut comparisons (96 overlapping-read cases; 24 inert controls). Destructive deletion control differs in maximal event count. Permanent hub growth proven; raw retained graph not qualified as active geometry. H1 restricted pass; no destructive-conflict completion. B2/B3 deferred. results/EX-GC31B1B/EX-GC31B1B_note.md. Review deferred by PI.
**EX-GC31B1C (ChatGPT):** two active edge replacement gadgets with immutable event ancestry, proof of fair causal completion and bounded degree, 240/240 exact normalized schedule/relabel checks; 224 endpoint-overlap cases. Full prior-edge turnover, max degree≤3 across tested seeds. Fixed refinement geometry; producer DAG is a forest with no joins, no branch interaction. H1 restricted pass; B2/B3 held. results/EX-GC31B1C/EX-GC31B1C_note.md. Review pending/deferred.
**EX-GC31B1D (ChatGPT):** adaptive local FIFO rendezvous, all 480 exact cut comparisons pass, concurrent adjacent events, genuine multi-parent cross-root joins, vertex/port-owner turnover, degree≤3/queue≤2. Marker alters active geometry away from source root in 3/5 seeds; no energy or metric claim. Permanent ≤3-port ancestry boundaries prevent unrestricted geometry. Fisher transformation and deterministic dataflow prior art documented; no novelty claim. results/EX-GC31B1D/EX-GC31B1D_note.md. B2/B3 held. Review deferred.
**EX-GC31B1E (ChatGPT):** native two-ended bond rendezvous, 144/144 exact schedule/relabel comparisons; full vertex and bond turnover, degree≤3, cross-root causal joins and inherited interface doubling; distributed mincuts double R2→R3. Ancestral adjacency invariant, cubic boundary exponent log2/log6; no new neighbor classes, no emergent geometry. results/EX-GC31B1E/EX-GC31B1E_note.md. B2/B3 held; review deferred.
**EX-GC31B1F (ChatGPT):** two-hop port-word partner discovery vs one-hop control, 240/240 exact cut comparisons, full turnover and max 4 incidences. New ancestral neighbors in 4/4 connected families; zero control/cross-component links; pair counts plateau R1–3. Distributed mincuts grow but prescribed routing/gadget geometry; no physical metric. results/EX-GC31B1F/EX-GC31B1F_note.md. B2/B3 held; review deferred.
**EX-GC32A (ChatGPT):** paper-only screen of seven constructions, including Sakharov/Jacobson/entanglement, causal sets/networks, Haramein and Heim. Exact stable Gaussian family has identical full ground state and total vacuum energy, microscopic exchange symmetry, but gapless cone ratios 1, 5/3, 3. Retarded impulse derivatives reconstruct the generator; fixed quadratic response is state-independent even when conserved excess energy increases. Geometry/clock supplied in this control; no candidate passes full gravity gates. [Note](results/EX-GC32A/EX-GC32A_note.md), [plan](spec/EX-GC32A_plan.md). Review/PI closure pending; negligible compute.
**EX-GC32B (ChatGPT):** 42/42 blind-to-adjacency graph/distance reconstructions, all 42 common-cone classifications correct. Same-state different speed and quadratic energy-independent response verified. 60 initial + 30 additional exact quantum nonlinear cases separate onsite vs edge effects; both gates pass. Same Hamiltonian and source density, total energy 1–3: first two response moments identical analytically, later response changes (at t=2, up to .0525 onsite/.2974 bond vs E=1, g=1). Density-mediated prescribed link response is not a total-energy source; full kernel sensitive to source dynamics. No geometry/gravity derived. <2 CPU s including rerun. [Note](results/EX-GC32B/EX-GC32B_note.md). Review/PI closure pending.
**EX-GC32C (ChatGPT):** NOT SUPPORTED for positive harmonic ZPF geometry selection. Exact monotonicity empty-graph minimum, all 32,768 labeled N6 graphs; 12/12 same-budget comparisons disconnected blocks < connected necklace < extended lattice. 48 bounded searches, all best graphs cluster >.1, 20–73% descending lattice swaps; no physical dynamics inferred. Excitation gradients independently verified, massive uniform mode positive E but zero link force; matched-E coherent preparations differ in force. Pure oscillator U0 insufficient, not a no-go for induced stiffness. 15.03 CPU s. [Note](results/EX-GC32C/EX-GC32C_note.md). Review/PI closure pending.
**EX-GC33A (ChatGPT):** paper-only screen complete; PI previously explored Graphity, no repeat. Positive-semidefinite link ground/free energy monotone and concave with fixed Hilbert/constraints. Exact C12 vs 2C6: both positive oscillator energy and half-filled fermion band favor disconnected comparator at all positive mu/weights. Minimal conversion resonance has energy split 2g with equal link occupancy; full-source operator identity derived for uniform lapse. Matrix geometry conditional scalar/spinor shared principal operator; cited specific IKKT one-loop EH action, not full metric/source derivation. [Note](results/EX-GC33A/EX-GC33A_note.md). Negligible compute, review/PI closure pending.
**EX-GC31A–EX-GC33A (11 studies: 31A, 31B1, 31B1B–F, 32A–C, 33A): PI CLOSED 2026-11-22** ("I sign off on GC31A–GC33A."). Claude review: review/EX-GC31A-33A_claude_review.md (math correct; series drifted from Model Gc). The "closure pending" wording in the per-study lines above predates this.
**GC33B NOT PURSUED** (PI decision 2026-11-22 after the Claude review). **Next: EX-GC34** (Claude leads; ChatGPT does separate recomputation/review): one shared set of link weights plus a site lapse used by every field (cubic + face-diagonal links encode a full spatial metric) on the GC26 knot. 34A common cone, 34B energy-source test (GC33A identity), 34C light bending / PPN gamma. Outline docs/EX-GC34_outline.md. PI go-ahead for 34A: "Go ahead with the 34A plan." (2026-11-22).
**EX-GC34A done (Claude, 2026-11-23), PI CLOSED 2026-11-23 (superseded by 34A-r):** registered centroid readout NOT MET (periodic-box wrap, my plan flaw); registered spectral recomputation shows one common cone for all fields (≤2e-5), control caught, bending ratio 2.009 (GR-like/lapse) and ~0 (Nordström). Options: 34A-r rerun with wrap-immune readout (~1.1 core-h, over cap) or go to 34B. results/EX-GC34A/EX-GC34A_note.md.
**EX-GC34A-r done (Claude, 2026-11-23), PI CLOSED 2026-11-23:** 68/69 registered cases pass (M1 54/54 one common cone; M2 14/15, F4 in GR-like bump just over tolerance). Light-bending ratio GR-like/lapse-only 2.005, Nordström ~0 on prescribed backgrounds. Next proposed: 34B (knot on shared links + source test), needs PI go-ahead. results/EX-GC34Ar/EX-GC34Ar_note.md; figure paper/figures/exgc34ar_cone.png.
**EX-GC34B done (Claude, 2026-11-23), PI CLOSED 2026-11-23:** lapse source = energy universally (incl. binding); static spatial source = lattice offset → 0; so light bending (PPN γ) is set by the gravity-field action, not the matter coupling. Proposed 34C: what in Model Gc selects that action (bootstrap/self-consistency), not yet approved. results/EX-GC34B/EX-GC34B_note.md.

**EX-GC34C done (Claude, 2026-11-23), PI CLOSED 2026-11-23:** γ = 1 (GR light bending) is selected by two rules together: "pure clock" (the lapse multiplies gravity's own energy, as for matter) and a single covariant gravity energy (√γR). Either alone leaves γ free. Suggested 34D (not approved): Regge-type lattice gravity energy, check γ ≈ 1 on the lattice. results/EX-GC34C/EX-GC34C_note.md.

**EX-GC34D done (Claude, 2026-11-23), PI CLOSED 2026-11-23:** Regge gravity energy (lapse × deficit angles, Kuhn triangulation) gives linearised GR's isotropic solution exactly on the lattice (γ = 1 at all k; lapse Newtonian with G = κ/8π). Registered packet-bending criteria P5/P6 FAIL (2.056/2.062 vs 2 ± 3%); post hoc isotropic-coordinate rerun gives 1.989/2.000, so a coordinate effect of the readout. Analytic: a moving source needs a shift-like variable. Suggested 34E (shift link variable, gravitomagnetism, wave speed) or a white-paper consolidation, not approved. results/EX-GC34D/EX-GC34D_note.md.

**EX-GC34E done (Claude, 2026-11-23), PI CLOSED 2026-11-23:** shift link variable + lattice DeWitt kinetic energy → linear GR on the lattice: gravitomagnetism matches GR (4e-4 at kh=0.05), two light-speed wave modes and nothing else (post hoc scan; registered wave-count method failed on a singular pencil), re-slicing symmetry exact along axes and O((kh)²) off. Suggested: white paper v0.5 consolidation of 34A–E with ChatGPT review, or 34F time-stepped waves. results/EX-GC34E/EX-GC34E_note.md.

**White paper v0.5 DRAFT (Claude, 2026-11-23, PI go-ahead "go ahead and do the white paper v0.5"):** `paper/Model_Gc_white_paper_v0.5.md` adds Section 9 (EX-GC27–33 summary) and Section 10 (EX-GC34A–E), Figure 10 (`src/whitepaper_figures_v05.py`). v0.4 unchanged. Not frozen: ChatGPT recomputation/review of EX-GC34 and of the draft pending; freeze needs PI approval.

**ChatGPT review of v0.5 and correction (Claude, 2026-11-23):** ChatGPT found a growing (imaginary-frequency) vacuum mode in EX-GC34E; Claude confirmed it post hoc (`src/exgc34e_review_unstable.py` → `results/EX-GC34E/review_unstable.json`): Im ω/k ≈ 1.42 (111, generic), 2.12 (110), 1.57 (near 100) at every kh from 0.025 to 0.2; lapse + lengthwise shift only (the re-slicing direction, ~99.98% coordinate change, curvature part ~0.06 (kh)²); the 100 axis undecided (degenerate pencil). The 34E "exactly two modes / linear GR" claims are withdrawn (note corrected); 34D note P5 corrected (packet 2.056 inside the window, rays 2.0615 outside). White paper v0.5 revision 2 adopts the review's verdict wording (still a DRAFT, not frozen). Proposed next: `spec/EX-GC34F_plan.md` (registered 2026-11-23, PI "run 34F") (after ChatGPT's plan review `review/EX-GC34F_draft_chatgpt_review.md`): diagnostics only, F1 spectrum without assumed symmetries, F2 constrained time-domain growth, F3 mechanism; a fix becomes a separate 34G. NOT approved. Unresolved: ChatGPT's 2.5e-4 wave/light speed gap vs our 1.4e-4.

**EX-GC34F done (Claude, 2026-11-23; PI go-ahead "run 34F"), PI CLOSED 2026-11-23:** the 34E growing re-slicing mode is real (constrained time evolution reproduces Im ω), persists as h → 0 at fixed k (Im ω/k direction-dependent 0.35–15.6; absent only on axes and (1,−1,0)-type diagonals, where re-slicing is exact). Cause: per-tetrahedron DeWitt kinetic term. Diagnostic: per-cube strain fit removes the growth and keeps the 34D static field and 34E gravitomagnetism, but adds a third slower mode; λ ≠ 1 also removes it (changes continuum). Several registered checks failed for method reasons (see note); post hoc argument-principle confirmation. Suggested next: 34G, a strain-rate construction that keeps exactly two modes and no growth (criteria in spec/EX-GC34F_plan.md); not approved. results/EX-GC34F/EX-GC34F_note.md.

**EX-GC34G done (Claude, 2026-11-23; PI go-ahead "Go ahead with 34G"), PI CLOSED 2026-11-23:** seven fixed kinetic-term candidates; none passes as registered. Cube-mean B is free of growth and of a third mode post hoc, but has strong k² anisotropy in one polarisation, static chains, and a non-dynamical edge combination (as do all fitted constructions) that broke the time integrator. Suggested next (not approved): either a registered follow-up on B (index-reduced integrator, chain gauge proof, anisotropy) or derive a kinetic term with exact lattice re-slicing. results/EX-GC34G/EX-GC34G_note.md.

**ChatGPT review of 34F/34G (2026-11-23, `review/EX-GC34F-34G_chatgpt_review.md`):** agrees with both verdicts. Caveats: 34F minimal model misses the generic growth rate by ~35% and local contours do not prove spectrum completeness; restrict conclusions to the tested discretisation. Priorities: (1) bounded analytical follow-up on B (identify the no-kinetic edge combination, derive its constraint, test whether the drift changes tidal curvature, constraint-reduced evolution, wave-speed error as a separate gate); (2) if that fails, derive a kinetic action from an exact lattice re-slicing identity; (3) keep the 34D/34F/34G distinctions in the white paper. Next step not chosen.

**EX-GC34H done (Claude, 2026-11-23; PI go-ahead "go ahead with the B run first."), PI CLOSED 2026-11-23:** constraint-reduced analysis of cube-mean B. B is not viable as registered (H2, H4, H5, H6 FAIL; H1, H3, C1 pass). The non-kinetic 7th-edge combination leaves a long Jordan chain in the zero sector: constraint-satisfying solutions grow polynomially and are mostly non-metric (body-diagonal mismatch); one polarisation up to 5.8% fast at kh=0.05 (∝ (kh)²). Suggested next (not approved): a 6-metric-component variable set (drop the body-diagonal edge from the dynamics) or a kinetic term derived from an exact lattice re-slicing identity. results/EX-GC34H/EX-GC34H_note.md.

**EX-GC34I done (Claude, 2026-11-23; PI "Go ahead with deriving a kinetic term from an exact lattice re-slicing identity"), PI CLOSED 2026-11-23:** requiring re-slicing to be an exact null vector of the kernel gives one linear identity K sigma_xi = C_en/l0; exact local (range-1) solutions exist with a forward-difference gradient, none with central. Registered If1 (closest to T): no growing root anywhere, exactly two isotropic low modes, cutoff third mode, chains = symmetries; fails G2 (1.1% slow, normalisation mu taken from T) and G5 (E4 2.3e-2). Post hoc mu = -1 passes G1-G3, G7; E4 1.25e-2, first order in kh; synchronous gauge gives a regular system with clean constraints. Suggested next (not approved): register mu = -1, cure the first-order shift response (staggered beta or wider central stencil), synchronous-gauge time domain with tidal observable. results/EX-GC34I/EX-GC34I_note.md.

**EX-GC34J done (Claude, 2026-11-23; PI "Go ahead with the suggested next step"), PI CLOSED 2026-11-23:** mu = -1 registered. Forward-gradient re-slicing-exact kinetic term (J1, range 1 closest to T) is stable, has exactly two isotropic waves at c (7.5e-4), cutoff 7th-edge mode, symmetric chains, clean synchronous-gauge time evolution (post hoc, constraint-filtered seeds). Fails only gravitomagnetism E4: first order (kh/4, 1.25e-2 at kh=0.05), identical for every exact solution, so forced by the forward gradient. Options (not approved): accept first order with finer lattice; staggered shift for both gauges; centred gradient at larger range or rational symbol. results/EX-GC34J/EX-GC34J_note.md.

**EX-GC34K done (Claude, 2026-11-23; PI "Maybe you should do all 3"), PI CLOSED 2026-11-23:** worst E4 case = source moving along a lattice axis, component (01|12), c1 = 1/4. The coefficient varies across exact kernels (91% spread, A-variable, narrowing 34J) but all 33 stay first order and none beats J1 on the max over cases. Gradient family: only forward/backward have local solutions (both first order); centred = averaged staggered shift is impossible for any bounded kernel (gradient vanishes at zone boundary where the lapse coupling does not, post hoc); intermediate and Kuhn gradients have no local solution (residual shrinks with range). J1 passes registered G4'' (clean seeds). Option (a): E4 <= 1e-2 at kh <= 0.03. Open (not approved): targeted minimisation of the max coefficient; non-local kernels; change the lapse coupling to admit a centred gradient. results/EX-GC34K/EX-GC34K_note.md.

**EX-GC34L done (Claude, 2026-11-23; PI "yes, go ahead. write the plan and run it."), PI CLOSED 2026-11-23:** per ChatGPT 34K review. Source derived from exact lattice continuity (link current) removes the kh/2 re-slicing mismatch exactly; gauge-corrected readout exactly invariant; E4 still first order (1.55e-2 at kh .05, ratio 2.00): L-FAIL. Linearised min-max of the first-order coefficient over all exact kernels (ranges 1, 2) only 0.310 -> 0.292, optimised kernels grow at kh = 3: C-FIXED. First order lives in the geometry (forward re-slicing gradient + vertex lapse coupling); 34K A2 spread was mostly readout/gauge dependence (inference). Remaining routes (not approved): accept first order (kh <= 0.03); change the potential side so a centred gradient is admissible. results/EX-GC34L/EX-GC34L_note.md.

**EX-GC34M done (Claude, 2026-11-23; PI "do number 1) Start with a scalar wave packet"), PI CLOSED 2026-11-23:** one discrete scalar action (P1 tet energy from edge lengths, lumped volume, vertex/edge lapse, forward-edge shift) varied wrt J1 lapse/shift/edges. Energy (re-slicing) identity exact; local momentum (displacement) identity first order, O(1) coefficient (M2 FAIL); total momentum and integrated stress fine (2nd order). Post hoc: moving momentum density to the vertex fixes momentum but breaks energy. One shift location cannot serve both J1 symmetries; suggests a staggered lattice where displacement also acts on link shifts (new vacuum action). Review step 2 should wait. results/EX-GC34M/EX-GC34M_note.md.

**EX-GC34P done (Claude, 2026-11-24; PI "go ahead with that test", ChatGPT 34O review), PI CLOSED 2026-11-23:** 34O scalar packet (normalisation fixed) sources Y1 linear gravity: constraints kept to RK4 error (34N stress breaks them 4-8%), tidal curvature matches linearised GR at 0.5-1.2% (h=.25), second order, incl. the velocity-dependent part. No feedback on scalar. ~2.1 core-h (just over cap). Next (not approved): feedback/self-consistent test, knot on Y1, origin of Y1 geometry, white paper. results/EX-GC34P/EX-GC34P_note.md.

**EX-GC34O done (Claude, 2026-11-24; PI "go ahead with 34O"), PI CLOSED 2026-11-23:** on Y1 the scalar's parallel normal stress in product form G_a(v)G_a(v-e_a) (action term 1/4 h_aa(v)G_a(v)G_a(v-e_a)) makes both conservation identities exact (momentum <= 2.7e-15, energy <= 4e-14) on packets and random fields; 34N's transverse and shear terms were already right. Next (not approved): review step 2, couple the packet to the Y1 vacuum and test whether its energy and momentum change gravity as linearised GR. results/EX-GC34O/EX-GC34O_note.md.

**EX-GC34N done (Claude, 2026-11-23; PI "redesign the shift couplings and vacuum action"), PI CLOSED 2026-11-23:** Y1 lattice staggered by index parity (n, phi, h_aa at sites; beta, u, momentum on links; h_ab on plaquettes); vacuum = linearised Einstein-Hilbert with staggered differences (k -> K = 2 sin(q/2)/h). Both gauges exact, range-1 local, no growth, exactly 2 waves at |K| (same cone as staggered scalar), gamma = 1, gravitomagnetism E4 1.04e-4 at kh .05 second order (J1 1.25e-2 first order), G4'' pass, packet energy exact, momentum ~2nd order (0.2-0.4% at h=.25) but registered M2 FAIL on coarse-step ratio. Spatial metric no longer Kuhn edge lengths. Next (not approved): finer-spacing M2 rerun; review step 2 (energy changes gravity) on Y1; knot on Y1. results/EX-GC34N/EX-GC34N_note.md.
*(Superseded by the PI decision above:)* **Next decision:** propose GC33B conditional matrix-geometry algebraic common-metric and complete stress variation test (scalar/spinor/gauge, metric-image rank and independent energy ledger). No GC33B authorized or run; GC31B2/B3 held. Prior PC allowance estimate 21–23 of 40 core-hours, not remeasured. Dated updates below are history.

## Set-up state (2026-11-18)
**Owner of the active work at set-up:** Claude.

## State
- The repository was set up on 2026-11-18. It holds the Gc1/Gc2 code, records and notes copied from ModelG (source
  commit edab7d47; see `PROVENANCE.md`), the model definition (`spec/MODEL_Gc.md`), and the protocol.
- Inherited results:
  - EX008–EX010 (exploratory).
  - MG195 (confirmatory, PI CLOSED 2026-11-18): P1 and P2 SUPPORTED, P3 UNRESOLVED, C2 PASS.

## Candidate next steps (none approved yet)
1. **Physics mapping note:** done, `docs/PHYSICS_MAPPING.md`. Suggested tests from it: weighted action-reaction with an unequal pair; drift bifurcation (inertia); telegraph G carrier (finite speed); G budget.
2. **Stochastic survival (EX-GC01):** noisy 3D continuum (tau-leap or Langevin), to find the minimum Ω at which a Gc
   particle survives. This is a prerequisite for a token version.
3. **Token / hypergraph Gc (2D, then 3D):** G tokens only hop and are destroyed only at W, so conservation is exact.
4. **Existence window** in (A, B, k, c), and the P3 retest with the short-range part separated (registered).
5. **Gc2 with built-in charge symmetry**, and the gravity/charge strength ratio.

## Update 2026-11-18: EX-GC01 done
Weighted action and reaction holds for particles with equal mobility μ; free fall is not universal when μ differs
(doubled k: 1.8 times faster). New design target: a universal μ across particle types (results/EX-GC01/EX-GC01_note.md).

## Update 2026-11-18: EX-GC02 done
Inertia appears when the core species W responds slowly (drag → 0 at τ_W ≈ 28; effective mass confirmed by a direct
time response), but the threshold is type-dependent, so free fall becomes non-universal. Revised design target: one
quantity setting both G consumption and inertia (results/EX-GC02/EX-GC02_note.md).

## Update 2026-11-18: EX-GC03 done
Push ∝ M (active = passive gravitational mass) but drag is set by the shared X/Y core, so free fall differs by
type (35%). Suggested next: a bound pair of identical particles, to see whether binding breaks universality
(results/EX-GC03/EX-GC03_note.md).

## Update 2026-11-18: EX-GC04 done
Identical particles bind at d = 7.145 and 10.11 (short-range Turing force), and a bound pair falls 6% (3%) more slowly
than a single particle. Cause: consumption grows with local G, so the upper partner pulls harder; the deficit fades as
(dM/dG)/(8πd). In Gc1 the same G-dependence that makes particles fall breaks the equivalence principle for composites
(results/EX-GC04/EX-GC04_note.md). Next idea: a variant where drift does not come from G-dependent consumption.

## Update 2026-11-18: EX-GC05 done (Stage A)
No rotating Gc1 particle: non-round core modes never oscillate; the X/Y chemistry oscillates only as a whole medium
(ether Hopf τ_Y = 8.64). Route 1 to a 720° particle needs a core with its own clock (e.g. a slow inhibitor of W);
routes 2 (orientation field) and 3 (hypergraph ribbons) need new ingredients too (results/EX-GC05/EX-GC05_note.md).

## Update 2026-11-18: EX-GC06 done
A slow local core inhibitor V ("Gc3") gives the particle a breathing clock (period ≈ 18–20) near the existence
edge (κ ≈ 0.9–1). Rotation does not come first: l = 1 oscillates only later, and l = 2, 3 never. Next would be a 3D
time simulation (does breathing saturate, does l = 1 become circling, any 2:1 lock) (results/EX-GC06/EX-GC06_note.md).

## Update 2026-11-18: EX-GC07 done
3D (spectral-in-angle) time runs of Gc3: breathing grows until the core collapses; no steady pulse seen, no circling.
Route 1 would need a saturating brake (another tuned ingredient); routes 2/3 get 720° from structure instead
(results/EX-GC07/EX-GC07_note.md). Solver src/exgc07_sim3d.py (8 ms/step) is reusable for other Gc variants.

## Update 2026-11-19: route 3 design note and EX-GC08 done
docs/ROUTE3_RIBBON_DESIGN.md: 720° from a ribbon-structured Gc hypergraph tethered by its G supply. EX-GC08 (Stage 1):
a local three-node rule (no edge may swing through the half-turn ambiguity) on cube-rotation frames reproduces the
belt trick (1 turn stuck, 2 turns untangle); pairwise rules fail. Next: Stage 2, Gc tokens on this ribbon graph
(results/EX-GC08/EX-GC08_note.md).

## Update 2026-11-19: EX-GC09 done (route 3, Stage 2)
Gc1 as tokens on a 3D cubic lattice (h = 0.5): the particle survives at Ω = 1000 (M −6%) and 300 (−20%) over T = 30,
fades at 100, dies at 30 (about 10⁴ core tokens needed). Far G deficit matches the lattice point-sink field (1–4%).
A splitting bias in the first token update was found and fixed (results/EX-GC09/EX-GC09_note.md). Next: Stage 3,
tethered particle (needs the PI's go-ahead).

## Update 2026-11-19: EX-GC10 done (route 3, Stage 3)
Particle tethered by its G supply (frames on the EX-GC09 lattice, no-cut lift rule, twisted edges pass G with
probability 1 − ε). One turn: twist stuck next to the core, index −1, M −7.9% (deterministic, 3 seeds); two turns
untangle, M within 0.6%. Controls (ε = 0, pairwise rule) as predicted. Token runs (Ω = 3000, partly on the PC) agree on
topology but the mass signal is noisy (shared drift artefact). Open: exchange sign, a real rotating particle.

Overall status summary as of 2026-11-19: docs/STATUS_2026-11-19.md.

## Update 2026-11-19: soliton route (EX-GC11 to EX-GC13)
EX-GC11: the frame-field tether gives no exchange sign (tethers lack identity). Richard chose the soliton route
(docs/SOLITON_PARTICLE_DESIGN.md). EX-GC12: a degree-1 knot of lattice frames is stable, degree conserved. EX-GC13:
loop classes (Pontryagin–Thom): 360° turn = 1, swap of two knots = 1, swap + turn = 0, degree 2 = 0; the loops run on
the lattice through allowed moves. Next (needs go-ahead): make the knot a Gc particle (consumes G, gravitates).

EX-GC14 (2026-11-19): knot as a Gc particle: the gated W core lives only on a knot; 1/r source confirmed; falling not resolved (lattice pinning at h = 0.5).

EX-GC15 (2026-11-20): at h = 0.25 the knot particle (λ = 1) falls at exactly the plain Gc rate (continuum within 2%); without the pull (λ = 0) the knot shrinks and the particle dies.

## Update 2026-11-20: EX-GC16 done (recoil rule)
Fixed-rate consumers that recoil from absorbed G tokens fall universally in the weak regime (5 types within 0.5% of
m₀E/3), attract each other with an inverse-square force (matches lattice prediction 1–11%), and break universality
by 2–7% when self-depletion is strong (m₀/Ω large). Note: results/EX-GC16/EX-GC16_note.md. Candidate next steps (not
approved): drag on a moving absorber (Le Sage objection); put the recoil rule on the knot particle.

## Update 2026-11-20: EX-GC17 done (drag vs inertia)
Head-on tokens give drag γ ≈ m₀/6 per unit, not inertia; the acceleration response is small and positive (anti-inertia,
lagging wake). Terminal fall speed (push / drag) is universal within ±2% for all types, even where the push alone
differs by 7%. So the token medium yields Aristotelian (overdamped) universal free fall. Newtonian inertia needs another
source (internal dynamics, or a medium with no rest frame as in Rueda–Haisch–Puthoff). Note: results/EX-GC17/EX-GC17_note.md.

## Update 2026-11-20: EX-GC18 done (recoil on the knot particle) — failed
Kicking only the knot tears it off the core; pushing W locally squeezes the core and does not move it; only a rigid
shift of the whole particle moves it (and then K 0.2 vs 0.4 differ by ~10%). A Gc particle is a chemistry pattern
with no way to carry momentum as a whole: the same gap as the missing inertia. Candidate next step (not approved):
give the Gc fields a velocity/momentum variable. Note: results/EX-GC18/EX-GC18_note.md.

## Update 2026-11-20: EX-GC19 done (flow/momentum field)
A shared flow field carried by the particle lets recoil move it as one body and accelerate it (K 0.2 vs 0.4 within 4%).
But the particle's own lagging G wake (G diffuses) acts as drag even with drag switched off (coasting time ≈ 46), and the
knot cannot keep up with the core at these speeds. Candidate next step (not approved): make G propagate as a wave
(finite speed) so the wake drag vanishes for uniform motion. Note: results/EX-GC19/EX-GC19_note.md.

## Update 2026-11-20: EX-GC20/21 done (wave-like G; faster knot)
Wave-like G reduces the kicked particle's speed loss only ~27% (0.023 → 0.017); a push-off control coasts exactly, so the
loss is a real backward push from absorbing G that is at rest (hypothesis: accretion drag in a medium with a rest frame).
Faster/stronger frame dynamics do not help the knot; carrying the knot with the flow keeps it within one lattice step.
Candidate next steps (not approved): coasting time vs particle strength (weak-field limit); Lorentz-invariant G field.

## Update 2026-11-20: EX-GC22 (ChatGPT, exploratory, complete locally)
Prescribed moving sink: steady wake drag tends to zero with K, while response per assigned inertia remains finite. Matched full Gc particle (no knot, 8 cases): K=0.001 retains 99.82–99.88% of initial flow over 20 units, but the W centroid slips and slows at h=0.5; physical coasting unresolved. D=80 reduces strong particle depletion by 61× but not late internal flow decay (0.0104 vs 0.00998). The post-hoc initial balance identifies a missing convective velocity derivative: changing consumption weights over a fixed radial V produces deceleration despite zero net recoil. Evaluating −v_boost ∂z V cancels 99.9% of initial strong-wave decay (not a dynamic test). **Recommended EX-GC23:** test advected V, conservative momentum balance, fine-grid centroid and flux response; separately confront the EX-GC17 bare absorption drag disabled here. Read `results/EX-GC22/EX-GC22_note.md` and `spec/EX-GC22_plan.md`. No claim of Newtonian inertia or Lorentz invariance. No study running; PI assigns next owner. ChatGPT local branch `ex-gc22-chatgpt` based on df7e73a; remote writes rejected (403). PI can import supplied bundle; no changes to remote main.

## Update 2026-11-20: EX-GC22 reviewed and merged (Claude)
EX-GC22's diagnosis is correct and reproduces (review/EX-GC22_claude_review.md). The EX-GC19/20 coasting loss was
mostly a flow-model defect (V not carried with the particle), so the "accretion drag" hypothesis is withdrawn (errata
added). **Proposed EX-GC23 (not yet approved):** add (V·∇)V, then a conservative momentum form. Measure the centroid
at h = 0.25 for the strong wave-like and a weak case, with flux trials and a larger box. Bare absorption drag
(EX-GC17) is to be handled separately afterwards.

## Update 2026-11-20: EX-GC23 (ChatGPT, exploratory)
Done: material transport fixes most earlier apparent coasting loss, fine-grid actual W follows flow within 2–3%. Strong K=.2 at h=.25, D=80, kick .05: W shift .96345 vs .84735 off; retained V 97.06% vs 73.91%; weak K=.01 on retains 99.81%, W shift .97476. Discrete momentum ledger closes and identifies cancellation of advected consumption reweighting. Flux K=.2/.4 actual centroid accelerations agree .36%; null stays stationary. Conservative alternative coasts similarly but M falls 4.87%, exact integrated cohesion zero; density still chemistry-assigned. Wider box ±9 agrees in W shift within .12% at t≈10, late residual decay still unresolved. All nine production runs complete, four independent settled fields, 1.616 core-h plus brief validation. See `results/EX-GC23/EX-GC23_note.md`, `spec/EX-GC23_plan.md`, `review/EX-GC23_response.md`. Branch `ex-gc23-chatgpt`; main unchanged, no study running. PI assigns Claude review and integration. **Candidate EX-GC24 (not approved):** reintroduce EX-GC17 head-on absorption drag in corrected fine-grid model and test explicit two-sided momentum exchange; inertia remains assigned, not derived. No claim of physical inertia or Lorentz invariance.

## Update 2026-11-20: EX-GC23 reviewed and merged (Claude)
EX-GC23 reproduces exactly (strong kick rerun bit-identical; own recomputation matches; review/EX-GC23_claude_review.md).
With (V·∇)V the particle really coasts: W shift 0.963 vs 0.847, flow kept 97% vs 74% (strong), 99.8% (weak); a ~2.5%
slip remains even when weak. Candidate EX-GC24 (not approved): reviewer suggests asking what sets inertia/drag (emergent
mass from medium momentum) rather than only reinstating the assigned drag; h = 0.125 slip check; orbits later.

## Update 2026-11-20: EX-GC24 candidate proposal (ChatGPT, at PI request)
Agrees with Claude: derive the inertia/drag ratio, rather than verify the assigned tau. Read `spec/EX-GC24_proposal.md` for the prescribed-sink impedance test (steady velocity, acceleration, harmonic response and parameter scaling), candidate weak-field added-mass expression, two-way momentum ledger and full-particle follow-up. Critical test: at fixed sink shape, the linear added-mass proxy scales as consumption squared, unlike linear G sourcing. Figures: `paper/figures/exgc24_particle_medium.png` (conceptual schematic, editable SVG beside it) and `paper/figures/ModelGc_topological_particle.png` (original topological-particle illustration). Neither is simulation data. No production approved or run; PI assigns execution. Orbits deferred.

## Update 2026-11-20: EX-GC24A (ChatGPT, exploratory)
A1–A4 completed after PI approval of Claude's reduced plan. Positive reactive inertia at weak damping, negative at strong damping; independent time integration agrees <.004%. S² scaling vs linear sourcing reproduced. With actual 1/G weighting and inherited head-on drag, reference m=.23124, b_total=5.34496, conditional m/b=.0433; weak .00217. A4 stops B under this closure; physical wave absorption law is not independently derived. A5 instantaneous grid diagnostic agrees with exact uniform boost but full .125 dynamic slip remains open (~14 core-hour estimated pair). See `results/EX-GC24/EX-GC24_note.md`, `spec/EX-GC24A_plan.md`, `src/exgc24_impedance.py`. Candidate next: design reciprocal Lorentz-compatible fields with conserved energy/momentum and stable boosted configurations or independently derive internal inertia. No orbits, no physical inertia claim. PI assigns next stage. No remote writes attempted.

## Update 2026-11-20: EX-GC25 (Claude, exploratory)
Short look at option (b), a medium with no rest frame (PI's choice after the EX-GC24A review). Argument: a conserved, absorbed carrier always defines a rest frame (uniform density = timelike four-vector), so the head-on drag can't be removed. A scalar field G with a wave equation keeps the 1/r field and 1/d² attraction with no rest frame; coupling it to energy (Nordström) gives charge = rest energy by von Laue's theorem, so free fall is universal and inertia is free. 1+1 D φ⁴ check: energy-coupled kinks with mass 8× apart fall together to 0.13%; non-energy coupling spreads them by Q/M; boosted kinks coast (loss ≤2e-4 over t=200) unless only ~2 lattice spacings wide. Costs: G no longer conserved/consumed (not Model Gc), particles must be conservative solitons (knot candidate), scalar gravity alone fails light bending/perihelion. See `results/EX-GC25/EX-GC25_note.md`. Candidate next (PI decides): EX-GC26, 3D relativistic knot with energy-sourced scalar G (pair 1/d², universal fall, coasting), or option (a) internal inertia.

## Update 2026-11-20: white paper drafts (Claude)
`paper/Model_Gc_white_paper_draft_v0.1.md`, `paper/Model_Gc_white_paper_draft_v0.2.md` (checkpoints) and the frozen `paper/Model_Gc_white_paper_v0.3.md` (through EX-GC25; second ChatGPT review applied). v0.2 which applies ChatGPT's review (endorsed by the PI): claim labels (demonstrated / imposed / interpretation), S vs M_i notation, drift-law wording for MG195, EX-GC16 as response per imposed inertia, classical-vs-quantum spin table, conditional EX-GC24A, gated seven-step direction. New figures `paper/figures/wp_loop_classes.png`, `wp_coasting_times.png` from `src/whitepaper_figures.py` (committed records only). Further review iterations expected.

## Update 2026-11-20: EX-GC26A (Claude, exploratory)
Plain Skyrme model on a lattice: a degree-1 Skyrmion is a stable, lossless particle (energy conserved to 3e-10 under a kick, degree holds, relaxes back to E₀ to 2e-7; no-knot control dissolves). 5 of 8 registered criteria met; not met: lattice degree at h = 0.2 (0.929), Derrick balance (0.911), lopsided start unconverged, momentum drift 1.3e-3 E₀ (lattice breaks translation invariance). Knot collapses through grids coarser than h ≈ 0.25. Compute ~2.4 core-h, over the 2 core-h cap (deviation logged). See `results/EX-GC26A/EX-GC26A_note.md`. Next (PI decision): 26B boosts at h ≤ 0.15, measuring lattice momentum leakage first; needs the PC.
## EX-GC26A2 (ChatGPT, exploratory; 2026-11-20)

Independent qualification complete, Claude review pending. See `results/EX-GC26A2/EX-GC26A2_note.md` and `spec/EX-GC26A2_plan.md`. Geometric degree exactly −1, lopsided energy agreement 8.73e-6; compact Fourier pinning tiny but raw pinning confounded by walls/interpolation. All six original-stencil short boosts pass 1e-3 E0 momentum gate; L=7 halves leakage; higher-order static candidate achieves Derrick ratio 1.00114 at h=.15 but not qualified dynamically. Compute under .37 core-h. Recommend gated exploratory 26B with larger box, mass tensor, small impulses and flux ledger; no precision inertia or gravity claim. PI decides. Original A scoring left to Claude.

## Update 2026-11-20: EX-GC26A2 reviewed (Claude)
EX-GC26A2 (ChatGPT) reproduces and checks out (review/EX-GC26A2_claude_review.md): geometric degree exactly −1, lopsided start converges, momentum leak mostly from the box. 26A scoring corrected to 4 of 8 fully met (C1 partial). Next (PI decision): 26B as a gated exploratory inertia check under ChatGPT's five conditions plus the review's additions (Lorentz-boosted initial data, stress-tensor momentum flux, optional pion-mass term to remove wall effects; target ~1%).
## EX-GC26B (ChatGPT, exploratory; 2026-11-20)

Complete; PI CLOSED 2026-11-26. Plan published before production d4c98a1. 15 dynamics cases, all endpoint geometric degrees −1, energy drift ≤3.66e-11. Fine h=.15,L=7 low-speed axial/diagonal mass errors +.734%/+.019%; axial .3/.5 misses 1% (+1.08%/+2.43%), transverse and high-speed flux residual gates fail. Larger box and half dt do not materially change .3 axial mass; spatial refinement halves errors. Centered static tensor max error .48%, matched stress caveat. Full precision gate not met; recommend low-speed exploratory continuation or targeted 26B2 stencil qualification before 26E. See `results/EX-GC26B/EX-GC26B_note.md`. 0.56 core-h production within cap.

## EX-GC26B review (Claude, 2026-11-21)

`review/EX-GC26B_claude_review.md`. All 15 runs recomputed exactly; axial/diagonal split reproduced with the original
26A NumPy solver (`review/EX-GC26B_check/`). P and H match Lorentz values within 0.85% at h = 0.15; the axial mass
excess is a low centroid speed: energy moves slower than P/H along the axes (0.971–0.993) but not diagonally (≥0.995),
consistent with the two-orientation stencil lacking single-axis mirror symmetry. Momentum loss at speed is lattice
non-conservation, not box flux. Recommend 26B2: eight-orientation stencil, speed/(P/H) and transverse-P gates, longer
coast, force-response test. PI CLOSED 2026-11-26.
## EX-GC26B2 (ChatGPT, exploratory; 2026-11-20)

Complete; PI CLOSED 2026-11-26. Plan published prospectively 897db1c. Eight-orientation mirror-symmetric solver eliminates transverse momentum (~2 million-fold improvement at .2/.3 axis), but primary v/(P/H) still fails .5%; exact semidiscrete bond energy-current expression explains old/new ratios to 1.7e-4. T5 speed loss .63% passes, global momentum drift 1.98e-3 fails. Linear weak-force pair M/E=.783, squared-deficit core control 1.056; both fail 1% force and impulse gates, strong tail/centroid coupling sensitivity. All ten dynamics endpoints degree −1, total energy drift ≤4.40e-11. .856 core-h. Precision 26E remains blocked; recommend current-consistent higher-order operator and force-response convergence tests. See `results/EX-GC26B2/EX-GC26B2_note.md`.

## EX-GC26B2 review (Claude, 2026-11-21)

`review/EX-GC26B2_claude_review.md`. All numbers reproduced. Claude's 26B mirror-symmetry explanation was only half
right (fixed transverse P, not the axial speed gap). ChatGPT's exact bond-current result confirmed independently on the
old stencil with the 26A code (`review/EX-GC26B2_check/`): J/P 0.9815 axis, 0.9988 diagonal vs fitted 0.9817/0.9988.
Post hoc: core-force P/v = 1.008 E0; the force miss is mainly dP/dt ≈ 0.96 fQ. Recommendation: stop stencil
qualification, record error budget (~1.5% h=.2, ~0.7% h=.15, ~h²), design 26E as a comparative test at matched w/h,
one short ramped core-force step first. PI CLOSED 2026-11-26.

## EX-GC26E1 (Claude, exploratory; 2026-11-21)

Done, PI CLOSED 2026-11-21. Plan 493cba7 before runs. B = 1 vs B = 2 torus in a ramped Nordström field (W = e^{2Φ} on the
quadratic term), h = 0.2, L = 7, same box, null-subtracted. Fall-rate ratio B2/B1 = 0.995–1.020 (baryon vs energy
centroid, two fit windows); lattice prediction (2E₂/E ratio) 1.013, continuum 1. Non-energy control: 1.345–1.445 vs
predicted 1.425. Binding 3.9%, so "binding does not gravitate" would give 1.041: suggestive (baryon centroid) but not
decisive (energy centroid 1.020). Momentum/impulse 0.964–0.969 for all runs (common lattice factor, = 26B3). 1.04 core-h.
Next: 26E2 at h = 0.15 (~3–4 core-h: cap increase or PC), then 26E3. See `results/EX-GC26E1/EX-GC26E1_note.md`.

## EX-GC26E2 (Claude, exploratory, PI's PC; 2026-11-21)

Done, PI CLOSED 2026-11-21. Plan 3ec3738; records 9e5eb3c (PC, MSVC build, PI option A), merged. h = 0.15: B2/B1 fall
ratio 0.988–1.009 (baryon/energy centroid, two windows), lattice prediction 1.006, continuum 1; no-binding alternative
2E1/E2 = 1.043 excluded by all readouts. Energy-centroid ratio moved 1.020 -> 1.009. Control 1.32–1.41 vs 1.39.
Momentum/impulse 0.982 (= 26B3). Candidate result: binding energy gravitates. 4.39 core-h of PC allowance (~29–31 left).
Next (PI decision): 26E3 (isospinning / breathing B = 1) or 26C (self-sourced field).

## EX-GC26C (Claude, exploratory; 2026-11-21)

Done, PI CLOSED 2026-11-21. Plan 2f8c24e. Knot + dynamical Nordström scalar (common action), h = 0.2, field box ±17,
G_N = 0.01 (and 0.02 for B = 1), T = 12. Field 1/r behind a light-speed front (rΦ flat to 2%); A/G_N = 0.995 × source,
→ 2E₂ as G_N → 0; B2/B1 active mass 1.936 (lattice 1.947, E ratio 1.922, no binding 2). Stable, energy 3e-12.
1.71 core-h cloud. Next (PI): 26D (boosted coupled knot) or 26F (two-body), better on the PC.

## EX-GC26D (Claude, exploratory, PI's PC; 2026-11-21)

Done, PI CLOSED 2026-11-21. Plan df01310 (amended ef866e8, 820ac4f before any PC run); PC records 6dec330. Boosted B = 1
knot with its own Nordström field (26C end state), u = 0.1, 0.3, T = 8, vs no-field controls. No field-induced drag:
coupled − control speed change −7e-5 / +5e-6. Energy ≤ 5e-10. The PC co-motion readout at d ≥ 3 was spoiled by the far
taper (my plan flaw); field-only check with exact far field (review/EX-GC26D_check) gives fore/aft 1.00–1.01 and
axial/side = 1/γ. PC 1.82 core-h (~27–29 of 40 left), cloud ~0.9 core-h. Next (PI): 26F (two-body) or 26E3.

## EX-GC26F (Claude, exploratory, PI's PC; 2026-11-21)

Done, PI CLOSED 2026-11-22. Plan d51adae; PC records 1b8ad8a. Two B = 1 knots (26C state, same orientation) attract
through their own Nordström field, d = 6 and 9.2, G_N = 0.01, T = 10, each minus a no-gravity control (direct Skyrme push
~1.5–2× gravity). Registered readout (baryon centroid) 0.901 and 0.960 of Newtonian prediction: SUPPORTED, narrowly;
all readouts 4–16% low (candidate: gravity × Skyrme-force cross term, untested). 1/d² favoured (2.20/2.38 vs 2.35;
1/d 1.53). Energy 2e-10. PC 6.42 core-h (~21–23 of 40 left). Next (PI): white-paper revision for 26A–F, deficit test, or 26E3.

## EX-GC26 series closed (PI, 2026-11-22: "I sign off: 26F and close the series")

26A–F are all closed. 26E3 (isospinning / breathing knot) was not run and would need a new go-ahead. Summary of the chain
(massless Skyrme reference model + Nordström scalar, not Model Gc): stable lossless knot (26A/A2), Lorentz-correct
moving knot with an O(h²) lattice P-vs-J gap (26B/B2/B3), binding energy falls with gravity (26E1/E2), knot sources a
1/r field with active mass = source (26C), coasts with its own field and no drag (26D), two knots attract ~1/d² at
0.90–0.96 of the Newtonian prediction (26F). Next (PI): white-paper revision for EX-GC26.

## EX-GC27A/B (ChatGPT, exploratory; 2026-11-22)

Complete, reviewed, PI CLOSED 2026-11-22. Registered plan `spec/EX-GC27AB_plan.md`. Static scale-sensitivity controls recover analytic Derrick scaling. Adding independent potential scale removes universality guarantee: virial-separated B2−B1 sensitivity −.002635/−.001775 at h=.2/.15; shrinking difference, composition dependence unresolved, no claim of new gravitational mechanism. All degree checks preserved, gradient validation 1.77e-9, 0.4012 core-h. `results/EX-GC27AB/EX-GC27AB_note.md` and summary. Recommend Candidate 4 / EX-GC28A normal-mode vacuum response as next high-value step; GC27C deferred until a motivated dynamical vacuum action exists. No quantum/dynamic tests run. PI decides.

## EX-GC26E3 (Claude, exploratory, cloud; 2026-11-22)

Done, PI CLOSED 2026-11-22. Plan 4c286fe, stage-2 registration c7e5563. Isospinning B = 1 knot (ω = 0.5, 0.7; kinetic
6–12%) in the 26E1 field. Registered baryon-centroid readout NOT SUPPORTED (falls 8%/13% faster than static, near the
"no extra inertia" alternative). Momentum ledger: weight = energy-coupling prediction to 0.1%, momentum/energy-centroid
velocity = total energy to 0.5%; the baryon centroid drifts away from the energy centroid when spinning. Untested
candidates: charge redistribution in the gradient, radiation from the sudden spin-up. 1.20 core-h cloud.

## EX-GC28A–28H and EX-GC29A (ChatGPT, exploratory; 2026-11-22) — PI CLOSED

Closed by the PI 2026-11-22 ("go ahead and close 28 and 29A and we will defer 29B"). Claude reviews (all runs
reproduced): `review/EX-GC27AB_28A_claude_review.md`, `review/EX-GC28B-G_claude_review.md`,
`review/EX-GC28H_29A_claude_review.md`; ChatGPT response `review/EX-GC27AB_28A_response.md`. Scalar-QFT, Skyrme and GR
reference models, not Model Gc. Lessons:
- 28A: radial-channel zero-point energy of the B = 1 Skyrmion is regulator-controlled once Born-subtracted; under a
  non-metric coupling its gravitational weight per unit energy (≈0.84) differs from the classical part (0.92).
- 28B: connected zero-point (Casimir) interaction of two defects in 3D falls as 1/d³ (force 1/d⁴): not Newtonian.
- 28C/28H: loop-induced collective scalars can be stiff and gapless only with tuning or a shift symmetry; the
  symmetric case forbids a monopole source, the dilaton-like case gives species-dependent sources (~3× spread).
- 28D: Sakharov induced Einstein–Hilbert term recovered with textbook coefficients; its size is set by the cutoff.
- 28E–G: with minimal metric coupling the whole dressed particle (binding included) gravitates via ADM mass, but a
  massless secondary scalar is a non-universal fifth force; screening reduces it only ~1/α and leaves species spread.
- 29A: particle size gives a generic elastic form factor (momentum-transfer roll-off), not a frequency cutoff, since
  forward scattering is never suppressed.
Open question carried forward: where a common metric coupling comes from (proposed EX-GC30). 29B deferred.

