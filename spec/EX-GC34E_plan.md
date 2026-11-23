# EX-GC34E — a shift link variable: does lapse × Regge plus shift give full linear GR on the lattice? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "you can close 34D and go on to 34E". Cloud, cap 2 core-hours
(expected < 0.2: Fourier-space linear algebra). Fifth stage of EX-GC34. Our construction, not LaViolette's Model G.

## Question

34D showed that a Regge gravity energy, multiplied by the lapse, gives GR's static weak field exactly on the lattice.
It also showed analytically that a moving source needs a shift-like variable. This study adds that variable and the
kinetic energy that goes with it, and asks two things:

1. **Gravitomagnetism.** Is the field of a moving mass right at first order in its velocity?
2. **Gravitational waves.** Does the lattice carry exactly two light modes moving at the speed of light, which is the
   speed every matter field on the same structure has (34A)? An extra light scalar mode would mean a non-GR theory.

## Construction (stated before any readout)

- **Variables per vertex:** lapse n, shift β (a 3-vector: the velocity of the lattice vertices relative to the
  clock-normal observers), and the 7 Kuhn edge-length changes δl of 34D.
- **Edge strain rate.** For an edge D (rest length l⁰ = |D|), the strain rate is the length change rate minus what
  the vertices' relative motion explains: Dᵀ K D = l⁰ (dδl/dt − D̂·(β_end − β_start)). That is the lattice form of
  2NK_ij = ∂_tγ_ij − D_iβ_j − D_jβ_i. In each tetrahedron the six edges fix the six K_ij uniquely.
- **Kinetic energy (pure clock: kinetic / N).** (1/2κ) Σ_tet V_tet (K_ijK_ij − K_ii²) / N, the lattice DeWitt form. At
  linear order N = 1 here.
- **Potential energy:** the 34D Regge term (1/κ) Σ N_e l_e ε_e, unchanged.
- **Matter:** dust with mass M_v and momentum P_v at vertices, coupled as −Σ N_v M_v + Σ β_v·P_v (the linear
  expansion of a point particle's −mN√(1 − (ẋ + β)²/N²)). The momentum-flux coupling ½h_ijT^{ij} is O(v²), and
  everything here is O(v).
- **Exact by construction:** vertex displacements u(t) with δl → δl + D̂·Δu and β → β + u̇ leave every term
  unchanged. They are fixed by a penalty, as in 34D. **Not built in:** the lattice version of re-slicing time
  (n → n + ξ̇⁰, β → β + ∇ξ⁰). In the continuum this is what removes the would-be third (scalar) mode. Whether the
  lattice keeps it is part of the question.

Smoke check allowed before registration: implementation only (displacement modes null, static limit = 34D kernel).

## Parts and registered criteria

**E1 — implementation checks**
- The vertex-displacement modes are null vectors of the full kernel M(ω, q) at sampled (ω, q): ‖MG‖/‖M‖ ≤ 1e-10.
- M(0, q) restricted to (n, δl) equals the 34D kernel to ≤ 1e-12.
- **Separate path:** the kinetic energy is evaluated directly in real space on a periodic 8³ lattice. For random
  complex single-frequency fields, K_ij is solved per tetrahedron and summed. This must match the Fourier quadratic
  form to ≤ 1e-10 relative.

**E2 — re-slicing (refoliation) symmetry, reported**
- Smallest singular value of the shift–shift block, and of the full kernel on the space of lattice-gradient shifts,
  relative to the largest, at several q. **Expectation:** this symmetry holds only approximately on the lattice,
  broken at O((qh)²). I'm not confident of this, and either answer is informative.

**E3 — gravitational waves**
- Vacuum modes come from the generalized eigenproblem after eliminating the shift through its own equation (the
  momentum constraint), fixing the displacement gauge, and keeping the lapse as a multiplier (the Hamiltonian
  constraint).
- **P1:** at kh = 0.05 and 0.1, along (1,0,0), (1,1,0) and (1,1,1), exactly two modes have ωh < 0.5.
- **P2:** their speed satisfies |ω/k − 1| ≤ 1e-3 at kh = 0.05, and it differs from the lattice light speed of a
  massless matter field on the same structure (34A isotropic split) by ≤ 1e-3.
- **P3 (reported):** any other finite modes, with their dispersion. Cutoff-scale lattice modes (ωh ≳ 1) are
  expected from the extra edge variable.

**E4 — gravitomagnetism (forced, first order in v)**
- Source: M̂ = 1, P̂ = M̂v, ω = k·v. Solve to first order in v around the 34D static solution. If the lattice
  operator is singular in some direction, use least squares and report the residual. Reconstruct h_00 = −2n,
  h_0i = β_i and h_ij (34D edge → metric map), then form the linearised Riemann tensor. Compare it with linearised GR
  for the same conserved source, solved independently.
- **P4:** for v parallel to k, perpendicular to k, and oblique (|v| = 1e-3), the O(v) Riemann components (including
  the gravitomagnetic R_0ijk) agree with GR to ≤ 1e-2 of their largest component at kh = 0.05, and the error
  shrinks with k.
- **P5:** for v ⊥ k (the case 34D's no-shift system got wrong), the solution is finite, the momentum constraint is
  satisfied (residual ≤ 1e-10), and P4 holds.

## Scope limits

Linear theory on a flat background only. The time step and any lattice time discretisation are not modelled
(continuous time). Matter is dust at O(v). This decides whether the shared structure with lapse, shift and Regge
links reproduces linear GR's vector and tensor sectors. It does not test non-linear dynamics.

## Files

`src/exgc34e.py` (kernel, checks, E2–E4; writes `results/EX-GC34E/analysis.json`), `src/exgc34e_check.py`
(real-space kinetic sum and independent GR solve; `check.json`). It reuses the 34D tetrahedron Hessians from
`src/exgc34d.py`.
