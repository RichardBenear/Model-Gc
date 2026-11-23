# EX-GC34M — sources from one discrete matter action: does a moving scalar packet satisfy J1's lattice conservation identities? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "Close 34K and 34L. I suggest that the next step is to do
number 1) Start with a scalar wave packet before returning to the knot. Derive its energy, momentum, and stress
sources by varying one discrete matter action with respect to J1's lapse, shift, and spatial variables. The first gate
is whether those sources satisfy the lattice conservation identities as the packet moves. GC34L used a prescribed
moving source, so this would test a missing link in the model."

The source is `review/EX-GC34J_followup_chatgpt_review.md`, recommendation 1. This plan and the `explore/LOG.md` entry
are pushed before any 34M computation. Cap: 2 core-hours on the cloud; I expect under 0.3. Our construction, not
LaViolette's Model G.

## Why this is the right first gate

J1's vacuum equations have two exact null identities: displacement (3 per vertex) and re-slicing (1 per vertex). So
the forced equations M x = f can be solved only if the source is orthogonal to both, i.e. W†f = 0. In real space these
are the lattice conservation laws.

**Energy, from re-slicing ξ** (n += ∂_t ξ, β_a += −(ξ(v+e_a) − ξ(v))/h):

  I_E(v) = −∂_t s_n(v) + Σ_a [s_β,a(v) − s_β,a(v − e_a)]/h = 0.

**Momentum, from displacement u** (β_a += ∂_t u_a, e_t += D̂_t·(u(v+d_t) − u(v))):

  I_P,a(v) = −∂_t s_β,a(v) + Σ_t D̂_t,a [s_e,t(v − d_t) − s_e,t(v)] = 0.

Here s = δS_m/δx are the matter sources. The force on the geometry is f = −κ s, consistent with 34E's
f = κ(M e_n − P e_β). The identities hold on shell if the matter action is invariant under the same lattice
transformations, together with a suitable transformation of the matter field.

34L used prescribed sources. Here the sources come from one action, and the packet obeys the field equation of that
same action.

## The matter action (one action; linear in x = (n, β, e), quadratic in φ)

- **Field.** A real scalar φ at the vertices, with mass m, on J1's Kuhn lattice (spacing h) in continuous time.
- **Building block.** Each Kuhn tetrahedron has a P1 (linear finite-element) gradient energy computed from its 6 edge
  lengths l:

    E_tet(l, φ) = ½ V(l) dφᵀ G(l)⁻¹ dφ,

  where G is the edge-vector Gram matrix built from squared lengths, V = √det G / 6, and dφ holds the vertex
  differences.
- **Edge weights.** On the flat Kuhn lattice this energy is exactly ½ Σ_edges w_e (Δ_e φ)², with w = h/6 per
  tetrahedron on axis edges and 0 on the face and body diagonals (checked before the plan). So the flat operator is
  the 7-point Laplacian.
- **Vertex volume.** V_v(l) = Σ over the tetrahedra at v of V_tet/4. It is h³ when flat.

**The action:**

  S_m = ∫dt { Σ_v V_v(l) [ ½(1 − n_v) φ̇_v² − ½(1 + n_v) m² φ_v² ]
            − Σ_tet [ E_tet(l, φ) + Σ_edges w_e ( ½ n̄_e (Δ_e φ)² + φ̄̇_e β_e Δ_e φ ) ] }

**Where each coupling sits:**
- **Lapse.** At the vertex in the kinetic and mass terms. Edge-averaged, n̄_e = (n_i + n_j)/2, in the gradient term.
- **Shift.** Projected on the edge: β_e = h β_a(v) for an axis edge from v along a. Only axis edges carry weight, so
  no path choice arises. The edge-averaged time derivative is φ̄̇_e.
- **Spatial variables.** Through l = l₀ + e. The stress source comes from ∂E_tet/∂l and ∂V_v/∂l.

**Why the energy identity should be exact.** At linear order in x, S_m is invariant under re-slicing combined with
φ → φ + ξ φ̇. The algebra uses Δ_e(ξφ̇) = ξ̄_e Δ_e φ̇ + φ̄̇_e Δ_e ξ, and each term cancels against the lapse or the shift
coupling.

**Why the momentum identity probably is not.** Displacement moves the mesh. A P1 field carried by the
piecewise-linear map is P1 on the moved mesh, but not on the original one. So no transformation of the vertex values
makes the action exactly invariant. I expect a violation that falls as the spacing squared. Exact discrete
translations still conserve the total lattice momentum.

**The packet.** It evolves by the flat field equation of S_m (x = 0):

  φ̈ = −m² φ + (1/h²) Σ_a [φ(v+e_a) − 2φ(v) + φ(v−e_a)],

solved exactly mode by mode by FFT on a periodic box. The initial condition is a positive-frequency Gaussian packet:
- mass m = 1;
- carrier |k₀| = 0.6 (group speed about 0.5 in the continuum);
- envelope width σ = 3;
- box L = 32 (physical units);
- direction axis (1,0,0) or generic (0.8, 0.5, 0.33);
- amplitude 1.

**Spacings:** h = 1, 0.5 and 0.25 (N = 32, 64, 128), with the same physical packet. Evaluated at t = 0, 6 and 12, by
which time the packet has moved about 6 units.

**Time derivatives of the sources** are computed exactly. The sources are quadratic in (φ, φ̇), so I polarise them with
(φ̇, φ̈), taking φ̈ from the field equation. There are no finite differences in time.

## Gates

- **M0, derivation check (registered).** On an N = 8 lattice with random φ, φ̇ and a random small x, the change in the
  Lagrangian, (L(ε x) − L(−ε x))/(2ε), must equal Σ s·x to a relative 1e-8. L is evaluated directly from tetrahedron
  geometry with perturbed lengths, a separate path from the analytic source formulas.
- **M1, energy identity (registered).** max_v |I_E| / max_v |∂_t s_n| ≤ 1e-10 at every spacing, time and direction.
- **M2, momentum identity (registered).** r_P = max_{v,a} |I_P| / max(max |∂_t s_β|, max |stress divergence|).
  - **M2-EXACT** if r_P ≤ 1e-10 everywhere.
  - **M2-SECOND** if not exact, but r_P falls by ≥ 3.5× per halving of h (h = 1 → 0.5 → 0.25) in both directions, at
    every time.
  - **M2-FAIL** otherwise.
- **M3, sources vs the continuum stress tensor (reported).** At each spacing I compare the integrated lattice
  sources with the continuum stress tensor of the same packet:
  - total energy;
  - total momentum;
  - the integrated stress tensor, mapped from the 7 edge sources to T^{ij} through the 34D edge-to-metric map.

  The continuum T_μν of a free scalar is evaluated from the band-limited (FFT) interpolant of the lattice field. I
  report the relative differences and their order.

**A packet "satisfies the lattice conservation identities"** if M1 and M2-EXACT hold. M1 with M2-SECOND means the
matter couples consistently to J1's re-slicing, and to its displacement symmetry up to second order.

## My expectations, stated before running

- M0 passes.
- M1 holds exactly. The re-slicing cancellation is algebraic.
- M2 is SECOND (not exact), for the mesh-motion reason above.
- M3 differences fall as h².

## Not in this study

There is no coupled solve and no curvature comparison. Review step 2 (energy changes gravity) and step 3 (geometry
redesign) depend on this gate.

## Files

`src/exgc34m.py` writes to `results/EX-GC34M/`:
- `m0` → `m0.json`;
- `run` → `identities.json` (M1, M2, M3).

It reuses `src/exgc34d.py` and `exgc34e.py` unchanged.
