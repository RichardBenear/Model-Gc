# EX-GC62 — What λ = 4 counts: a written analysis of the shift strength, with the entanglement angle

2026-11-28 | EXPLORATORY, paper only (no runs) | Claude | PI request 2026-11-28 17:28 UTC ("it should be a note in the repo. You can do the D study after that"). Prompted by Cao, Carroll and Michalakis, *Space from Hilbert Space: Recovering Geometry from Bulk Entanglement* (arXiv:1606.08444), which the PI and Codex found in the literature. Our construction; not LaViolette's Model G and not Chain Theory.

## Summary

**Main point (section 1b).** λ is not independent: once γ is fixed, the absence of a preferred frame forces λ = 2(1 + γ). This is the PPN condition α₁ = 0, and it is what 52 and 53 measured (λ* ≈ 3.91). So the factor to count is 1 + γ, which 37B does in 3D. The remaining 2 is what local boost invariance demands. A counting rule that would supply it: the shift is the forward-minus-backward Doppler asymmetry of the **signal** (light) count.

- **λ = 4 is two factors of 2, and they are different kinds of thing.**
  - The first 2 is light against the clock: 1 + γ.
  - The second 2 is the shift against light. It is the same in every dimension.
  - In D space dimensions linearised GR gives exactly **λ = 2(1 + γ) = 2(D − 1)/(D − 2)**, so 4 only in three dimensions.
  - Our counting so far (37B) produced the first factor. Nothing has produced the second.
- **Measured in units of the spatial ruler, which is what an entanglement or area count measures directly, every coupling is a whole number:**
  - the clock is D − 2;
  - light is D − 1;
  - the shift is 2(D − 1).

  In 3D these are 1, 1, 2 and 4. On a cubic lattice, 2(D − 1) is the number of plaquettes that share a given link, and D − 1 is the number of coordinate planes containing it.

  So one candidate answer to "how many things are counted" is: **momentum registers on every plaquette around a link (4 in 3D), and light registers once per plane through the link (2 in 3D).** This is a hypothesis drawn from the pattern of integers, not a derivation.
- **37B's rule gets the right numbers only because D = 3.** Clocks are one-site events and signals are two-site coincidences, so light/clock = 2 in any dimension. GR needs 3/2 in D = 4 and has no Newtonian clock coupling at all in D = 2. 54's Doppler-squared count gives λ = 4 in any dimension, which GR rules out for D ≠ 3. (Correction, same day: that failure is inherited from the two-site signal weight, i.e. from 37B, not from 54's Doppler step; see section 1b.)
- **What this means for the substrate.** In D = 3 the clock and the ruler couple equally (γ = 1), so a clock-first substrate (event rates, ours) and a ruler-first substrate (entanglement and area counts, Cao et al.'s) cannot be told apart there. They separate in D = 2 and D = 4. That is the D study the PI asked for next.
- **The entanglement route to the second 2.** Cao et al. derive only the energy equation, G₀₀ = 8πG T₀₀. The shift equation is the momentum constraint, G₀ᵢ = 8πG T₀ᵢ, and they assume it by covariance.
  - Jacobson-style, applying the entanglement first law for every boosted local frame brings in the cross term 2kᵗkⁱT₀ᵢ. That is the same symmetric double count (g₀ᵢ dt dxⁱ + gᵢ₀ dxⁱ dt) as the factor 2 above.
  - Our 38B already showed the boost (Bisognano–Wichmann) form of the lattice entanglement Hamiltonian, so this step is within reach.

## 1. The exact statement from linearised GR in D + 1 dimensions

Use the trace-reversed perturbation h̄_μν = h_μν − ½η_μν h, in Lorenz gauge, with □h̄_μν = −16πG T_μν. Take a slowly moving source: T₀₀ = ρ and T₀ᵢ = −jᵢ, with stresses negligible. For d = D + 1 spacetime dimensions, h_μν = h̄_μν − η_μν h̄/(d − 2), and h̄ = −h̄₀₀. This gives:

| Quantity | Relation to h̄₀₀ | Field law |
|---|---|---|
| clock (lapse) Φ, with g₀₀ = −(1 + 2Φ) | h₀₀ = h̄₀₀ (D − 2)/(D − 1), Φ = −h₀₀/2 | ∇²Φ = 8πG (D − 2)/(D − 1) ρ |
| ruler ψ, with gᵢⱼ = (1 − 2ψ)δᵢⱼ | hᵢⱼ = h̄₀₀/(D − 1) δᵢⱼ | ∇²ψ = 8πG/(D − 1) ρ |
| light, the potential a ray feels (Φ + ψ) | h₀₀ + hₓₓ = h̄₀₀ | ∇²(Φ + ψ) = 8πG ρ |
| shift β = g₀ᵢ | h₀ᵢ = h̄₀ᵢ, with no trace part | ∇²βᵢ = 16πG jᵢ (up to sign convention) |

Hence:
- γ = ψ/Φ = 1/(D − 2);
- 1 + γ = (D − 1)/(D − 2);
- **λ ≡ (shift coupling)/(clock coupling) = 2(D − 1)/(D − 2) = 2(1 + γ).**

| D (space) | γ | 1 + γ (light/clock) | λ (shift/clock) | shift/light |
|---|---|---|---|---|
| 2 | — (clock uncoupled) | — | — | 2 |
| 3 | 1 | 2 | **4** | 2 |
| 4 | 1/2 | 3/2 | 3 | 2 |
| 5 | 1/3 | 4/3 | 8/3 | 2 |

D = 2 is 2 + 1 gravity: a mass makes a cone (ψ), with no Newtonian pull (Φ). That is a useful check that the table is right.

Two details matter:
- The light coupling, 8πG, does not depend on D. Trace reversal splits h̄₀₀ between the clock slot (share D − 2) and each spatial slot (share 1), and a light ray adds them back up.
- The shift is not split, because g₀ᵢ has no trace part.

So **shift/light = 2 in every dimension.** The 2 is the shift entering the line element twice, as g₀ᵢ dt dxⁱ + gᵢ₀ dxⁱ dt, whereas a potential enters a rate (a square root of the interval) at half weight.

## 1b. λ is not a third independent number: given γ, the absence of a preferred frame fixes it

(This was added after the first draft. It changes what has to be explained.)

- **Boost argument.** Take a static source's field, h₀₀ = −2Φ and hₓₓ = −2γΦ, and boost it to velocity v. To first order the boosted field has

  |h₀ₓ| = v |h₀₀ + hₓₓ| = 2(1 + γ) |Φ| v.

  So a moving source's shift is the Lorentz image of the light potential (h₀₀ + hₓₓ), and its strength relative to the clock is 2(1 + γ). The factor 2 is the boost mixing t and x in both slots of the symmetric metric.
- **PPN form.** In the parametrised post-Newtonian metric, the transverse part of g₀ᵢ carries the coefficient 2(1 + γ) + α₁/2, where α₁ measures preferred-frame effects. **α₁ = 0 holds exactly when λ = 2(1 + γ).**
- **Our own data say the same.** 52 found a preferred frame with no shift (α₁ ≈ −8). 53 removed it with λ = 4, and the combined 52 + 53 data put the required value at λ* ≈ 3.91. That is the measured form of this identity, with γ = 1 from 37B's counting.

**Consequence.** The "how many things are counted" question for λ splits cleanly:
1. **1 + γ** is a genuine counting question about the clock and signal counts. 37B answers it in 3D. The dimension dependence sits entirely here.
2. **The factor 2** is not a further count to discover. It is what local boost invariance demands of whatever carries the shift. The mechanism question is then: **what substrate process makes a moving source's field the boosted image of its static light potential?**

**Counting reading of that process.** For a source moving at v, each endpoint's tick rate, seen along a link, is Doppler-shifted by (1 + v) or (1 − v). If the shift is the forward-minus-backward asymmetry of the **signal** count (the light potential, weight 1 + γ), each endpoint of the coincidence is Doppler-shifted, so each direction's count changes by (1 ∓ v)² ≈ 1 ∓ 2v. The one-way delays per unit length are then −(1 + γ)Φ(1 ∓ 2v). Their difference, 4(1 + γ)Φv, equals 2h₀ₓ, so the shift is 2(1 + γ)Φv. That gives λ = 2(1 + γ) with nothing tuned. (Correction: the first draft used a single Doppler factor per direction, which gives only (1 + γ) = 2 in 3D. The square, one factor per endpoint, is what supplies the boost's factor 2.)

**Correction (same day).** The first draft said 54 anchored the shift to the clock count. That was wrong. 54's D2 rule is this rule: each direction of a link counts its own two-site (signal) coincidences, Doppler-shifted at both endpoints, (1 ∓ v)², giving 2 × 2 = 4. In 3D the two are identical. 54's miss was in the operator that stands for "count × velocity": the lattice energy current ½{Ĥ, v̂} gives strength 3.85 in uniform N (against 53's momentum law at 4), and 3.37 inside the packets' own clock-rate wells, because its hops carry N weights. 54 itself named the open question: count per local tick, which removes the well factor.

So the specific prediction is: **shift = directed asymmetry of the signal count, normalised by the same coincidence weights that set light's speed.**

## 2. Counting in ruler units

The ruler ψ is the natural unit for an entanglement substrate. In Cao et al. the area of a region is half the sum of the mutual information on the edges it cuts (their eq. 9 and 16). An energy perturbation changes that count and so produces spatial scalar curvature, R = 2G₀₀ (their eq. 60). In other words, entanglement counts the ruler directly.

In ruler units the couplings are:

| Thing | Coupling / ruler | D = 2 | D = 3 | D = 4 | Lattice count that matches |
|---|---|---|---|---|---|
| ruler (spatial link length) | 1 | 1 | 1 | 1 | — |
| clock (site rate) | D − 2 | 0 | 1 | 2 | directions normal to a plane (cubes per plaquette ÷ 2) |
| light along a link | D − 1 | 1 | 2 | 3 | coordinate planes containing the link |
| shift on a link, per unit momentum | 2(D − 1) | 2 | 4 | 6 | plaquettes containing the link (both sides of each plane) |

The last column is a reading, not a theorem. It is supported in one place: on the lattice, the operator that the momentum constraint applies to the shift is the lattice curl-curl, (Cᵀ C β)_link. That is a sum over exactly the 2(D − 1) plaquettes that contain the link; the 58 curl operator C has this structure. So the shift field law genuinely "asks" every plaquette around a link. Whether the source does too, and with what weight, is the open question.

How it maps onto what our substrate can count:
1. **One-site events (clock ticks),** at rate N. They are affected by energy, through Φ.
2. **Symmetric two-site coincidences (signals),** at rate N_x N_y. 37B: these give the light/clock factor 2, which equals (D − 1)/(D − 2) only in D = 3.
3. **Directed coincidences,** x→y against y→x. A moving source alters each endpoint by about (1 ± v), so the forward-minus-backward difference is the only count that is odd in v. That makes it the natural place for the shift (53–58 put the shift on links as a directed transfer).
4. **Plaquette circulations,** the curl around a face. They are affected by the curl of the momentum (58's spin events live here).
5. **Cut edges of the entanglement graph,** the area or ruler count (Cao et al.).

The question "how many things are counted" becomes concrete: in 3D, does a unit of momentum register on 4 things per link (the plaquettes) while a unit of energy registers on 2 per link (the planes), with a clock reading 1? And in D dimensions, does the pattern become 2(D − 1), D − 1 and D − 2?

## 3. Our earlier attempts, read this way

- **37B: clock = 1 site, signal = 2 sites.** This gives light/clock = 2 in any D. It is right in 3D only. Read in ruler units, two sites in 3D equals two planes through a link, so 37B cannot tell which of the two is being counted.
- **54: Doppler-squared,** (1 ∓ v)² → 2 × 2 = 4. One 2 is light's two-site weight (1 + γ in 3D). The other is the Doppler square, one factor per endpoint, which is the boost's factor 2 and holds in any D. On our substrate the two-site weight is 2 in every D, so the rule gives 4 everywhere; the D ≠ 3 failure is 37B's, not the Doppler step's. Its numerical miss in 3D (3.85 in uniform N, 3.37 inside the wells) comes from the N weights in the energy-current operator. (Corrected; the first draft said it counted the wrong things.)
- **53–61:** λ = 4 is put in by hand. All the consistency results (no preferred frame, spin sourcing, local transport, the probe reading free-space GR) hold with that value, so the target is well defined.

## 4. The entanglement angle and the second 2

Cao et al. obtain G₀₀ = (2πκ/α)⟨T₀₀⟩ from the first law δS = δ⟨K⟩ for regions at rest, then covariantise by assumption (their eq. 62). The shift lives in the part they assume.

Jacobson (1995; 2016 entanglement equilibrium) obtains the full equation by demanding the first law for every local boost direction. For a null or boost vector k, the source term T_μν k^μ k^ν contains the cross term 2k⁰kⁱT₀ᵢ, so the momentum is counted twice, once as T₀ᵢ and once as Tᵢ₀. This is the same factor 2 as in section 1.

On our lattice:
- 38B showed that the half-space entanglement Hamiltonian is the boost generator (Bisognano–Wichmann, 0.8% to 0.02%).
- So the measurable statement is: put a moving source next to an entanglement cut, and the change in the modular Hamiltonian contains the momentum with weight 2 relative to the energy's light-coupling term.

That is a lattice computation, not an assumption, and it would test the "second 2" from the entanglement side.

## 5. What the D study should test (to be planned next)

These predictions hold for any rule that is to be GR in D + 1 dimensions:

| Ratio | D = 2 | D = 3 | D = 4 |
|---|---|---|---|
| light/clock = 1 + γ | — (no clock coupling) | 2 | 3/2 |
| shift/light | 2 | 2 | 2 |
| shift/clock = λ | — | 4 | 3 |
| shift/ruler | 2 | 4 | 6 |

Candidate rules and what each predicts:

| Rule | light/clock | shift/light |
|---|---|---|
| 37B clock-first (site events, two-site coincidences) | 2 in any D (wrong in D ≠ 3) | — |
| 54 Doppler-squared | — | gives λ = 2 × (two-site weight) = 4 in any D; wrong in D ≠ 3 only through 37B's weight |
| Ruler-first or plaquette counting (this note's reading) | (D − 1)/(D − 2) | 2 |

The ruler-first rule's numbers come straight from the integer counts in section 2. The study would build the lattice counts in D = 2, 3 and 4 (cheap, because the lattices are small and these are static source solves) and read the four ratios. The rule passes only if it hits all of them without D-specific tuning.

## Caveats

- Section 1 is standard linearised GR (Lorenz gauge, slow source). The shift law in transverse gauge has the same coefficient. The signs of T₀ᵢ and β are convention-dependent; only magnitudes are used.
- The integer pattern in section 2 could be a coincidence of small numbers. That is why the test must use D ≠ 3, where the counts differ.
- Cao et al. work with spatial geometry only, with no dynamics and no Lorentz invariance. Their energy is modular energy, related to T₀₀ only at a conformal fixed point. Nothing here relies on their specific Φ(I) distance function.
- λ = 4 relative to our κ assumes, as 51–61 do, that our N is the lapse with Φ = N − 1 and that β = g₀ᵢ.
