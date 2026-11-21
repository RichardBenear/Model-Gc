# Handoff — current state

## Current state (2026-11-20)
**No study is running.** The research through EX-GC25 is written up in the frozen white paper
`paper/Model_Gc_white_paper_v0.3.md` (reviewed twice by ChatGPT, endorsed by the PI). **Next step (PI decision):**
the seven-step gated sequence in its Section 10, starting with a conservative 3D topological particle (EX-GC26), with
acceptance criteria fixed in the log entry before running; or the alternative of internal inertia within Model Gc.
The white paper is not to be revised again until EX-GC26 gives results. The dated updates below are the history.

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

Complete, review/PI closure pending. Plan published before production d4c98a1. 15 dynamics cases, all endpoint geometric degrees −1, energy drift ≤3.66e-11. Fine h=.15,L=7 low-speed axial/diagonal mass errors +.734%/+.019%; axial .3/.5 misses 1% (+1.08%/+2.43%), transverse and high-speed flux residual gates fail. Larger box and half dt do not materially change .3 axial mass; spatial refinement halves errors. Centered static tensor max error .48%, matched stress caveat. Full precision gate not met; recommend low-speed exploratory continuation or targeted 26B2 stencil qualification before 26E. See `results/EX-GC26B/EX-GC26B_note.md`. 0.56 core-h production within cap.

## EX-GC26B review (Claude, 2026-11-21)

`review/EX-GC26B_claude_review.md`. All 15 runs recomputed exactly; axial/diagonal split reproduced with the original
26A NumPy solver (`review/EX-GC26B_check/`). P and H match Lorentz values within 0.85% at h = 0.15; the axial mass
excess is a low centroid speed: energy moves slower than P/H along the axes (0.971–0.993) but not diagonally (≥0.995),
consistent with the two-orientation stencil lacking single-axis mirror symmetry. Momentum loss at speed is lattice
non-conservation, not box flux. Recommend 26B2: eight-orientation stencil, speed/(P/H) and transverse-P gates, longer
coast, force-response test. PI closure pending.
## EX-GC26B2 (ChatGPT, exploratory; 2026-11-20)

Complete, Claude review and PI closure pending. Plan published prospectively 897db1c. Eight-orientation mirror-symmetric solver eliminates transverse momentum (~2 million-fold improvement at .2/.3 axis), but primary v/(P/H) still fails .5%; exact semidiscrete bond energy-current expression explains old/new ratios to 1.7e-4. T5 speed loss .63% passes, global momentum drift 1.98e-3 fails. Linear weak-force pair M/E=.783, squared-deficit core control 1.056; both fail 1% force and impulse gates, strong tail/centroid coupling sensitivity. All ten dynamics endpoints degree −1, total energy drift ≤4.40e-11. .856 core-h. Precision 26E remains blocked; recommend current-consistent higher-order operator and force-response convergence tests. See `results/EX-GC26B2/EX-GC26B2_note.md`.

## EX-GC26B2 review (Claude, 2026-11-21)

`review/EX-GC26B2_claude_review.md`. All numbers reproduced. Claude's 26B mirror-symmetry explanation was only half
right (fixed transverse P, not the axial speed gap). ChatGPT's exact bond-current result confirmed independently on the
old stencil with the 26A code (`review/EX-GC26B2_check/`): J/P 0.9815 axis, 0.9988 diagonal vs fitted 0.9817/0.9988.
Post hoc: core-force P/v = 1.008 E0; the force miss is mainly dP/dt ≈ 0.96 fQ. Recommendation: stop stencil
qualification, record error budget (~1.5% h=.2, ~0.7% h=.15, ~h²), design 26E as a comparative test at matched w/h,
one short ramped core-force step first. PI decision and closure pending.

## EX-GC26E1 (Claude, exploratory; 2026-11-21)

Done, PI closure pending. Plan 493cba7 before runs. B = 1 vs B = 2 torus in a ramped Nordström field (W = e^{2Φ} on the
quadratic term), h = 0.2, L = 7, same box, null-subtracted. Fall-rate ratio B2/B1 = 0.995–1.020 (baryon vs energy
centroid, two fit windows); lattice prediction (2E₂/E ratio) 1.013, continuum 1. Non-energy control: 1.345–1.445 vs
predicted 1.425. Binding 3.9%, so "binding does not gravitate" would give 1.041: suggestive (baryon centroid) but not
decisive (energy centroid 1.020). Momentum/impulse 0.964–0.969 for all runs (common lattice factor, = 26B3). 1.04 core-h.
Next: 26E2 at h = 0.15 (~3–4 core-h: cap increase or PC), then 26E3. See `results/EX-GC26E1/EX-GC26E1_note.md`.
