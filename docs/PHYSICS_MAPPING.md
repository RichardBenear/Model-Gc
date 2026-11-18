# What Model Gc's ingredients would correspond to in physics (note, 2026-11-18)

Author: Claude, at the PI's request (Richard Benear, 2026-11-18: "Let's do the mapping to physics first"). This is
analysis and interpretation only: there are no new runs. Each item is marked **established** (it follows from the
equations or from a run already recorded), **testable** (a prediction we can check in Gc), or **speculative**. Model Gc
is our construction, not LaViolette's.

## Summary
In its steady, far-field form Model Gc is **Newtonian gravity written as a conserved sink flow**. G∞ − G plays the role
of the gravitational potential, the consumption rate M plays the role of mass, and Gauss's law is the conservation of
G. Three things do not map onto known physics as the model stands:
1. **No inertia.** Particles move with a velocity proportional to the field gradient, not an acceleration.
2. **Diffusive propagation.** Changes in the field spread like heat, not at a finite speed.
3. **A one-sided budget.** Matter only consumes G, so the cosmic supply must be topped up or run down.

Each of these is a definite question for the next studies, not a reason to stop. The same structure applied to charge
(Gc2) has no budget problem, because + and − particles make and consume C equally.

## 1. The field: Poisson's equation (established)
Outside the cores, the steady G equation is D_G∇²G = k·G·W, and W is zero there, so ∇²G = 0. Define a potential
Φ = −α(G∞ − G), with α > 0 a unit conversion. Then

  ∇²Φ = (α/D_G)·s(x),  with s = k·G·W the local consumption density,

which has the form of Newton's ∇²Φ = 4πG_N·ρ, with ρ ↔ s and 4πG_N ↔ α/D_G. A particle's "gravitational mass" is its
total consumption M = ∫s dV. The far field is −M/(4πr), and the flux of G into any sphere around the particle equals M
whatever the sphere's size: that is Gauss's law, and here it is just G conservation. EX008 confirms this to 1e-6.

**What it means.** A conserved carrier with sinks at matter and supply from far away *is* the standard way to write
the Newtonian field. That is why Model G (whose ether makes and destroys G everywhere) is screened, and Gc is not.

## 2. Motion: velocity, not acceleration (established; the main gap)
A Gc particle in a G gradient moves at v = μ∇G. This is overdamped motion, like a bead in syrup: there is no
momentum. Two consequences:

- **Free fall is universal (established at the 1–2% level).** EX010 found the same μ, within 1–2%, for +, − and
  neutral particles. If μ is universal, every particle moves at the same speed in a given field, the overdamped
  analogue of the equivalence principle. A particle's own speed in a pair depends only on its partner's M:
  v₁ = μM₂/(4πd²).
- **Action and reaction balance in a weighted sense (testable).** From the formula above,
  M₁v₁ = μM₁M₂/(4πd²) = M₂v₂. So the "M-weighted centre" of a pair stays fixed, exactly like a centre of mass in
  Newton's third law, if M plays the role of inertial mass. This is a sharp prediction for a pair of *unequal*
  particles (different k or W strength). It holds only in the far-field linear regime, and only if μ is the same for
  both.
- **Inertia (speculative route).** Newtonian motion needs an acceleration, so speed is remembered. Reaction–diffusion
  spots can acquire this. Near a "drift bifurcation", a spot's velocity becomes a slow variable with its own relaxation
  (Krischer & Mikhailov 1994; Bode et al. 2002). Spot dynamics then becomes second order: velocity persists, and a
  force changes it. Finding that regime in Gc, for example by slowing the inhibitor (W or Y), is the most direct path
  from "falls at a speed" to "falls with an acceleration". Without it, orbits are impossible: two particles simply
  spiral together.

## 3. Propagation: diffusion, not waves (established; second gap)
G obeys a diffusion equation, so a change at a source is felt at distance r after a time of about r²/D_G, and with no
sharp front. Gravity changes travel at c (GW170817 constrains the speed of gravitational waves to about 1e-15 of c).
- **A minimal fix (speculative):** give G a finite response time τ. This gives the telegraph equation
  τG_tt + G_t = D_G∇²G − s. Signals then travel at c = √(D_G/τ), and the steady (Poisson) field is unchanged.
- **What that would map to:** a scalar field obeying a wave equation with matter as its source is **Nordström's scalar
  gravity** (1913). It has the correct Newtonian limit but fails two classic tests: it gives no bending of light, and
  the wrong sign for Mercury's perihelion shift. So a single conserved scalar can at best reproduce Newtonian gravity.
  Matching general relativity would need the carrier to shape the geometry itself (a tensor field, or in the
  hypergraph picture, G consumption changing the graph's connectivity). That is a long-term question, not a near-term
  test.

## 4. The global G budget (established; third gap, with a cosmological reading)
Every particle consumes G continuously, so a universe of N particles consumes NM per unit time. Two options:
- **A cosmic source** that replenishes G uniformly, at a rate that exactly matches total consumption on large scales.
  A uniform source would reintroduce a weak ether term. If it is spatially uniform and independent of local G, it does
  not screen: it adds a constant to s, which cancels against the mean consumption, like the "Jeans swindle" that
  Newtonian cosmology uses. This is LaViolette's continuous-creation picture in a different place: the ether supplies
  G, but not locally in response to the particle.
- **G∞ running down.** M is proportional to G near the core, and the pair force scales with M (and with μ, which also
  depends on the background), so falling G∞ would make gravity weaken over cosmic time. That is a varying-G cosmology
  in Dirac's sense. Lunar laser ranging limits |Ġ_N/G_N| to about 1e-13 per year, so any running-down would have to be
  extremely slow compared with the particle's own consumption time.

**Charge has no such problem (established for Gc2's structure).** + cores make C and − cores consume it. A neutral
universe (equal numbers) has zero net C production, so C∞ stays fixed with no cosmic source. This gives a natural
asymmetry: charge can be exactly balanced, gravity cannot, because there is no negative mass. In Gc the reason is that
G is only ever consumed.

## 5. Strength: gravity can be weak without being short-range (established)
EX009 showed that the force scales as k² while the field stays 1/r down to k = 0.001. In Model G, a weak G coupling
also meant short range (screening length about 1/√σ, EX005–EX007). In Gc, strength and range are independent, so a
gravity 10⁻³⁶ times weaker than the charge force is not ruled out by structure; it only needs k_G ≪ the charge
coupling. Whether one chemistry gives that ratio naturally is the open Gc2 question.

## 6. Pictures this resembles, and where it differs
- **Ether-sink gravity** (Riemann, Bjerknes, Yarkovsky; LaViolette's own G-on consumption). Gc is the version where the
  sink picture is actually consistent with the medium's equations: there is no hidden ether rebalancing, unlike Model
  G (see `docs/lineage/`).
- **Le Sage gravity** (matter absorbs a flux of particles). Gc is effectively Le Sage with a *diffusing* medium instead
  of a ballistic one. That removes Le Sage's drag problem (there is no preferred direction in a diffusive bath) but
  keeps its heating problem: the consumed G becomes X, which the particle must shed (X → Ω). In Gc that energy input
  is part of what sustains the particle, which is the dissipative-structure view.
- **The "river" model of general relativity** (space flowing into masses, Painlevé–Gullstrand coordinates). The analogy
  is only loose. In the river model the inflow speed grows as 1/√r and its flux is not conserved, whereas Gc's G flux
  is conserved, like the Newtonian field lines. Gc maps onto the Newtonian potential, not onto the river velocity.

## 7. Tests this note suggests (none run yet)
1. **Weighted action and reaction** with an unequal pair: does M₁v₁ = M₂v₂ hold (§2)? This would be a short
   extension of the MG195 pair solver to two different particles.
2. **Search for inertia**: is there a drift bifurcation in Gc (a slow inhibitor or W time scale) where particles move
   on their own and respond to a gradient by accelerating (§2)?
3. **Finite propagation speed**: a telegraph G carrier. Does the particle survive, and does the far field still reach
   1/r with a retarded front (§3)?
4. **Budget**: G∞ falling slowly in time, or a uniform cosmic source. How do M, μ and the pair force follow (§4)?

References: Krischer & Mikhailov, Phys. Rev. Lett. 73, 3165 (1994); Bode, Liehr, Schenk & Purwins, Physica D 161, 45
(2002); Nordström, Ann. Phys. 42, 533 (1913); Abbott et al. (LIGO/Virgo), ApJL 848, L13 (2017); Williams, Turyshev &
Boggs, Phys. Rev. Lett. 93, 261101 (2004).
