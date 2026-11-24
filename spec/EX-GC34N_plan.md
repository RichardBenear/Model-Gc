# EX-GC34N — a staggered lattice where lapse, shift, metric and matter all sit where both gauge symmetries act (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "close 34M and then redesign the shift couplings and vacuum
action".

This plan and the `explore/LOG.md` entry are pushed before any 34N computation. Cap: 2 core-hours on the cloud; I
expect under 0.3. Our construction, not LaViolette's Model G.

## The design problem (from 34J–34M)

J1 keeps the shift β and the displacement u at the vertices, but its re-slicing acts through a forward difference,
which lands half a cell away.

- 34M showed that a matter action whose energy is exactly conserved must put its momentum density on the links.
- J1's displacement symmetry wants that momentum at the vertex.
- Averaging between the two placements gives a centred gradient. 34K showed that this has zeros at the zone
  boundary, where J1's lapse coupling does not vanish.

So no shift placement works on J1's lattice. The fix has to move the shift, the displacement and the metric
components together.

## The new layout ("Y1": staggered by index parity)

Each field sits at the site offset given by the parity of its spatial indices, as in the Yee/MAC staggering used in
electromagnetism and fluids:

| field | continuum meaning | lattice position |
|---|---|---|
| lapse n, scalar matter φ, re-slicing parameter ξ | h₀₀ = −2n, scalars | sites v |
| shift β_a, displacement u_a, matter momentum density | h₀a, vectors | links v + e_a/2 |
| normal metric h_aa | h_aa | sites v |
| shear metric h_ab (a ≠ b) | h_ab | plaquette centres v + (e_a + e_b)/2 |

Every derivative in a gauge transformation is then a centred one-cell difference between neighbouring positions:
- n += ∂_t ξ;
- β_a += −Δ_a ξ/h + ∂_t u_a;
- h_aa += 2 Δ_a u_a / h;
- h_ab += (Δ_a u_b + Δ_b u_a)/h.

Each is exact and local. In Fourier space, with amplitudes referred to each field's own position, every Δ/h becomes
i K_a with K_a = 2 sin(q_a/2)/h.

**The vacuum action.** It is the linearised Einstein–Hilbert (Fierz–Pauli) action in continuous time, with each
spatial derivative replaced by the corresponding staggered difference. In Fourier space it is the continuum quadratic
form with k_a → K_a.

- **Why exact gauge invariance survives.** The continuum form is gauge-invariant as an algebraic identity in the
  wavevector, so the lattice action is exactly gauge-invariant under all four transformations.
- **Why it is local.** Each term is a full index contraction, so the half-cell offsets pair up into a finite-range
  stencil.
- **What changes from J1.** The 7 Kuhn edge lengths are replaced by 6 staggered metric components. The non-metric
  7th-edge mode therefore disappears.
- **What stays the same.**
  - The lapse still multiplies the linearised scalar curvature: 34C's "pure clock" structure, and a single covariant
    gravity energy.
  - Only q = 0 has K = 0 in the zone, so there are no doubler modes.

**The scalar matter action** (Part E). It is 34M's action moved to the staggered layout:
- φ at sites; gradients G_a = Δ_a φ / h on links.
- Kinetic and mass terms at sites: ½(1 − n + ½ tr h) φ̇² and −½(1 + n + ½ tr h) m² φ².
- On links:
  - gradient term −½(1 + n̄ + ½ tr h̄ − h̄_aa) G_a² (bars are averages from the two adjacent sites);
  - shift term −φ̄̇ β_a G_a.
- On plaquettes: shear term + h_ab Ḡ_a Ḡ_b (bars are averages of the two bordering links).

The energy identity should still be exact, by 34M's algebra. Momentum density and displacement are now on the same
links, so I expect the momentum identity to hold at second order.

## Gates

**A. Construction**
- **I0.** For 20 random (q, ω), the 4 gauge vectors are null vectors of the 10 × 10 kernel M(q, ω) to ≤ 1e-12
  relative, and M is Hermitian (≤ 1e-14).
- **G7, locality.** M(q + 2π e_a) = M(q) to rounding. The real-space stencil from an N = 16 FFT vanishes beyond range
  2 (≤ 1e-13 relative).

**B. Spectrum (synchronous gauge n = β = 0; 6 metric components; companion matrix)**
- **G1.** No root with |Im ω| > 1e-8 |K| on the 34F/34G direction set at kh = 0.0125, 0.025, 0.05, 0.1, 0.2, plus
  axis, (1,1,0) and (1,1,1) at kh = 0.5, 1, 2, 3 and the zone boundary.
- **G2.** Exactly two wave pairs, ω = ±|K|, and their speed within 1e-3 of 34A's light speed at kh ≤ 0.05. Reported:
  the wave speed equals the massless staggered scalar's speed (same K).
- **G3.** At ω = 0 the static kernel has nullity exactly 4 (the static gauge) at every tested q ≠ 0: no extra
  static zero modes.

**C. Static field and gravitomagnetism**
- **G5-static.** A static point mass gives γ = 1, with the transverse projector built on K, at kh = 0.05, 0.5 and 2
  (≤ 1e-10 from 1).
- **G5-E4.** The 34K/34L moving dust source, 9 cases, kh = 0.0125, 0.025, 0.05 and 0.1:
  - density at the site, momentum on the links;
  - momentum p_a = ρ v_a k_a/K_a, which satisfies lattice continuity exactly;
  - Riemann tensor read with K_μ = (−ω, K);
  - compared with linearised GR at k.

  Gate: E4 ≤ 1e-2 at kh = 0.05 and order ratio E4(0.05)/E4(0.025) ≥ 3. Reported: the Riemann change under each lattice
  gauge (exactly invariant expected).

**D. Time domain**
- **G4″.** As in 34K: synchronous gauge, generic and (1,1,0) directions at kh = 0.05, 40 periods, two steps.
  - Seeds are constraint-compatible from the outset: all constraint-satisfying non-zero modes, plus the
    constraint-satisfying part of the zero cluster, each at 1e-8.
  - Pass: constraint residual ≤ 1e-9, wave tidal amplitude change ≤ 1e-6, each seed group ≤ 1e-6 of the wave, and
    both steps agree.

**E. Matter (34M's gates on the new lattice)**
- Same packet as 34M: m = 1, |k₀| = 0.6, σ = 3, L = 32; axis and generic; h = 1, 0.5, 0.25; t = 0, 6, 12.
- **E-M0.** A direct finite-difference check of the Lagrangian against Σ s·x, Richardson extrapolated, ≤ 1e-8.
- **E-M1.** The energy (re-slicing) identity, ≤ 1e-10 relative.
- **E-M2.** The momentum (displacement) identity:
  - **EXACT** if ≤ 1e-10;
  - **SECOND** if it falls by ≥ 3.5× per halving of h, in both directions at every time;
  - **FAIL** otherwise.

**Overall.** Y1 passes if I0, G7, G1, G2, G3, G5 (static and E4), G4″, E-M0 and E-M1 hold, and E-M2 is EXACT or
SECOND.

## My expectations, stated before running

Everything passes:
- E4 is second order, around 1e-4 at kh = 0.05;
- E-M2 is SECOND.

The risk I see is in E: the shear coupling averages the gradients to the plaquettes, which might leave a first-order
momentum mismatch in the off-axis directions.

## Not in this study

- The coupled matter–gravity solve (review step 2).
- The knot.
- Any change to the white paper.

## Files

`src/exgc34n.py` writes to `results/EX-GC34N/`:
- `build` → `build.json` (A);
- `spec` → `spec.json` (B);
- `static` → `static.json` (C);
- `sync` → `sync.json` (D);
- `matter` → `matter.json` (E).

It reuses `exgc34e.gr_solution` and `riemann` for the continuum comparison, and the 34F/34G direction set.
