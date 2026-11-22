# EX-GC34 outline — one shared link structure for all matter (DRAFT, not registered)

2026-11-22 | Claude | proposal for PI decision | EXPLORATORY track. Not a frozen plan: a plan per stage goes into
`spec/` and `explore/LOG.md` before any run.

## Question
Can gravity in Model Gc come from **one** set of link variables that every field must use to move, so that
(1) all particles and waves share one light cone by construction, (2) energy sources those link variables through one
action, and (3) the result bends light, which the EX-GC26 scalar (Nordström) field cannot?

## Why this, and why now
* GC30–33 showed a common cone never comes for free when each sector has its own coefficients. The working
  exception is when one structure builds every kinetic term (GC33A's matrix-model observation).
* EX-GC26 already had half of this: one action, W = e^{2Φ} multiplies every matter term, so fall is universal and
  binding gravitates. But Φ is a single scalar per site. A scalar conformal factor cannot bend light, and the
  coupling was put in by hand.
* Model G link: in SQK the G concentration sets local clock and transport rates. Here a site "lapse" plays that role,
  and the new ingredient is direction-dependent link weights (the spatial metric).

## Lattice ingredient
On a cubic lattice the 3 axis links give only a diagonal metric. Adding the 6 face-diagonal directions gives
9 link directions per site, enough to encode a full symmetric 3×3 spatial metric (6 components) plus redundancy to
test isotropy. Every field's gradient energy is Σ_links w_ℓ (Δ_ℓ field)²; its time term carries the site lapse N.
No species-specific coefficient anywhere.

Fields: the GC26 Skyrme knot (matter, B = 1 and B = 2), a massless test wave χ ("light"), and a massive test field ψ.

## Stages (each gated on the previous)
**34A — common cone is structural (cheap, cloud, < 0.5 core-h).** Fixed w and N backgrounds: uniform, anisotropic,
and a smooth bump. Measure the dispersion and fronts of χ, ψ and small pion waves on the knot field. Pass: one
characteristic cone for all three in every background, within lattice error. Negative control: the GC31A
two-coefficient case, which must fail. Expected to pass. Its job is to validate the code and remove "tuned speeds"
as an objection.

**34B — what energy sources the links (cloud or PC, about 1 core-h).** Give w and N the simplest local wave/stiffness
action (no Einstein structure built in). The source is then the matter's own link energies, which is the lattice
stress tensor. Apply GC33A's source test, total source = κ × total energy, to a knot at rest, a B = 2 torus (binding),
a boosted knot (kinetic) and a χ wave packet (radiation). Prediction to register: rest, binding and boost pass
through the lapse, and radiation sources differently (trace-free). Report it whichever way it comes out.

**34C — light bending, the headline (PC, about 2 to 4 core-h).** A knot with its own dynamical w and N, and a χ
packet passing at impact parameter b. Measure the deflection against two references: Nordström (0) and GR
(4GM/bc², twice Newtonian). Equivalently, read off PPN γ from the spatial part of the knot's field. Possible outcomes:
γ = 0 (scalar again), γ = 1 (GR-like), or something else. If γ = 1 appears only when the w action is given gauge
structure (linearised Einstein form), record that as the finding rather than tune toward it.

**Later, not in this proposal:** a token/rewrite version of the shared-hopping rule to restore the causal-invariance
requirement, and GC29B.

## Stop rules
Stop at 34A if a common cone needs per-species adjustment. Stop at 34B if the w action is unstable or the knot is
torn apart. Don't tune the w action to hit γ = 1.

## Roles and compute (suggested)
Claude leads: it reuses the GC26 knot code (`exgc26c.py`, `exgc26e1_stencil.cpp`) and is Claude's proposal.
ChatGPT does the separate recomputation and review of each stage. Total under about 6 core-hours, staged, and each
stage inside its own cap.
