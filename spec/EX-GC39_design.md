# EX-GC39 design note — does matter set the lapse, or does the lapse shape matter?

2026-11-26 | DESIGN NOTE (no runs) | Claude | Requested by Richard (2026-11-26: "I'm not sure that matter sets the lapse. It could be that the lapse sets how matter emerges. Go ahead and do a design note that tests the idea you have."). Nothing here is approved to run. Each study below needs its own plan, LOG entry and go-ahead. This is our construction, not LaViolette's Model G.

## The two hypotheses

Results so far:
- 37B: a local firing rate N (the lapse) gives the weak-field metric with γ = 1.
- 37C: N must be a relativistic wave field, not a diffusing carrier.
- 38A: the vacuum's entanglement gives the infrared part of gravity's stiffness.
- 38B: the vacuum looks thermal with respect to a boost clock at every cut, and a lapse moves that temperature the Tolman way, at about 60% of the predicted size.

None of these says which way the causation runs.

- **H-A, matter sets the lapse.** N (with the spatial metric) is fixed at each moment by the matter around it, through a constraint. Jacobson's version: the geometry adjusts so that the entanglement of every small region stays at its vacuum value ("entanglement equilibrium"). N has no life of its own. Remove the matter and N relaxes to 1. There are no free scalar gravitational waves.
- **H-B, the lapse shapes matter.** N is the primary field, with its own dynamics. Matter is a pattern that forms in, and is held together by, the local clock-rate structure. A particle's existence, size and mass are set by N. Then N must be able to carry energy and waves with no matter present.

These are not exclusive. GR has both a constraint (the G₀₀ equation fixes the spatial curvature from matter) and matter living in a given geometry. What separates the hypotheses is **which of the two is fundamental**. That gives three measurable signatures:
1. Is the geometry *determined* by the matter's entanglement (A), or free (B)?
2. Does matter built from the substrate's own events keep its proper size and clock rate in a lapse? This is the equivalence principle, which B needs if N is to shape matter without being detectable locally.
3. Does N radiate on its own? B with a wave N predicts scalar gravitational radiation. GR, and so A, predicts none.

There is already one data point. **37A-S2 found that Y1's vacuum gives the lapse no stiffness**: the q² part of Π_nn is −9e-6 against Π_nT = 0.043. The vacuum alone does not make N a dynamical field. That leans toward A, but only on Y1 and only for the induced part.

## 39A — entanglement equilibrium (tests H-A)

**Question.** For a static lump of matter on the lattice, is there a spatial metric that keeps the entanglement of every small block at its vacuum value? If so, does that metric obey Newton's equation with a definite G?

**Setup.** The 3D lattice scalar of 38A, with the metric coupling of 37A-S2: kinetic 1/√γ, gradient √γγ^{ab}, mass √γ m².
- Conformal perturbation h_aa = 2ψ(x), i.e. γ = (1 + 2ψ)δ.
- Gaussian states, with blocks B of 2³ to 6³ sites at many positions.
- Entropies in double precision. The entanglement Hamiltonian K_B is computed in arbitrary precision, as in 38B.

**Steps.**
1. **Geometry response.** Compute δS_B/δψ(x), the change in each block's vacuum entanglement when the metric is perturbed. Expect δS_B ≈ 2η·δ(Area) from the area law, with η the 38A coefficient. This checks that the lattice entanglement "sees" area.
2. **Matter response.** Use two kinds of matter, because they behave differently:
   - **(a) A coherent lump** (a classical wave packet, or a static displaced field). This leaves every block's entropy exactly unchanged (δS_B = 0), because the displacement is a product of local unitaries. But it raises δ⟨K_B⟩.
   - **(b) A massive defect** (a local δm²), which changes the state itself, so δS_B ≠ 0.
   - Jacobson's argument uses δ⟨K⟩ as the matter term. Whether the lattice should use δS or δ⟨K⟩ is a real question, and (a) and (b) separate the two.
3. **Equilibrium solve.** Find ψ(x) such that δS_B[ψ] + (matter term)_B = 0 for all blocks, by least squares over all block positions and sizes. Read off the profile, and G_eq from the 1/r far field: ψ → G_eq M/r, with M = ∫δ⟨ℋ⟩.

**Predictions (to register in the plan).**
- **P1:** a solution exists with residual small compared with the matter term, and with the attractive sign (ψ > 0 near the mass, as GR's isotropic γ = 1 + 2GM/r).
- **P2:** G_eq = 1/(4η), where η is the 38A area coefficient, within 10%. This is the internal consistency of the equilibrium idea.
- **P3 (descriptive):** G_eq against 37A-S2's induced static coupling. 38A suggests they differ by 2–5× at the cutoff.
- **P4 (descriptive):** the coherent lump against the defect. Does the coherent lump, which carries δS = 0, produce any metric response at all? If the lattice "uses" δS, coherent matter would not gravitate. That would be a direct failure of A, since a classical field certainly carries energy.

**What it means.** If P1 and P2 pass, the spatial geometry is fixed by matter through the vacuum's entanglement, which is H-A. With γ = 1 (37B), the lapse follows. If no consistent ψ exists (large residual, or block-size dependent), entanglement equilibrium does not determine geometry on this lattice, and A loses its mechanism.

**Cost.** Small: blocks of at most 216 sites, a few hundred positions. Estimated under 1 core-h. It reuses `exgc38a.py` (entropy) and `exgc38b.py` (arbitrary-precision K).

## 39B — matter emerging in a lapse (tests H-B)

**B1. Equivalence-principle scaling of emergent particles.** If N shapes how matter forms, a particle in a region of lower N must still have the same *proper* size and *proper* frequency. Otherwise local physics would differ from place to place, and particles would fall non-universally (the EX-GC01–04 failures). 37B's rule fixes how each kind of event scales. One-site events (clocks, local reactions) run at rate N. Two-site coincidences (hops, propagation) run at rate N². GR's isotropic metric needs coordinate size ∝ N and coordinate period ∝ 1/N.
- **Reaction–diffusion particle (Gc1, `exgc09_lattice.py`).** Reactions at rate N, diffusion by hops at rate N². Then the coordinate size goes as √(D/k) ∝ √N, so the proper size goes as N^{-1/2}. Expected to **fail** the equivalence principle.
- **Wave-built particle (Skyrmion, GC26).** On-site potential terms at rate N, gradient terms carried by coincidences. Then the coordinate size goes as N²/N = N, and the proper size is constant. Expected to **pass**.
- **Test.** Run both in uniform N = 0.8, 0.9, 1.0 and 1.1. Measure coordinate size, internal frequency and rest energy. Then put each in a gentle N gradient and measure the fall rate against the universal one.
- **Predictions:**
  - **P1:** the Skyrmion keeps its proper size and frequency within 1%;
  - **P2:** Gc1's proper size scales as N^{-1/2} within 10%, and it falls non-universally by the corresponding amount.
- **Meaning.** The scaling argument is nearly analytic, so this is a check, not a discovery. Its value is that it turns "the lapse shapes matter" into a hard constraint: **matter must be made of the same coincidence events as light**, as 37D already hinted. Diffusive matter cannot be shaped by N without breaking the equivalence principle.
- **Cost.** Under 1 core-h, using existing code.

**B2. Can a lapse well hold its own matter?** This is the strongest version of B: a particle is a self-consistent lapse well plus the matter mode trapped in it, with nothing else binding it.
- **Setup.** 37C's wave-equation N, sourced by the energy of a massive lattice field. The field lives in H_N, as in 38B.
- **Steps.**
  - Find the bound modes of a given lapse well, i.e. the threshold depth for a mode to appear.
  - Then iterate to self-consistency: occupy the mode, compute its source, re-solve N.
- **Questions.**
  - Does a stable, localized, self-bound lump exist?
  - Is there a discrete mass spectrum, which would mean N selects particle masses? Or a continuous family?
- **Expectation, recorded.** The continuum analogue is the boson star, a continuous family whose size scales as 1/(G m²). With G at the cutoff scale (38A), the lumps would be cutoff-sized: "particles from the lapse" at the Planck scale, not at the masses we see.
- **Prediction (P3):** the existence of the lump and the shape of its mass–size relation. Scored against the Schrödinger–Newton scaling, M·R ∝ 1/(Gm²), within 10%. A deviation would be the interesting result: the lattice selecting masses.
- **Cost.** About 1 core-h (3D eigenproblem in a box).

## 39C — does N radiate on its own? (separates A from B observationally)

If N is a free wave field (H-B), an oscillating or orbiting pair of sources radiates scalar waves. GR (H-A for the scalar sector) radiates only through the tensor modes, at the quadrupole rate. Binary pulsars confirm the GR rate to about 0.1%.
- **Setup.** 37C's comoving wave-N solver, made time-dependent. Two sources on a circular orbit in a box with absorbing walls.
- **Readout.** The radiated power against orbital frequency and separation, compared with the scalar quadrupole (and any dipole term) of a Nordström-type field.
- **Prediction.** Under B, nonzero scalar flux, scaling as the scalar quadrupole formula. Under A, none. A scalar flux of the computed size would rule out "N alone is a free field" observationally, unless N's waves are somehow screened.
- **Caveat.** Our N is a scalar, so it cannot produce GR's tensor waves either way (noted 2026-11-26). This test asks only whether N *adds* radiation that is not observed.
- **Cost.** 1–2 core-h.

## Order and decision points

1. **39B1 first.** It is cheapest and directly tests your version of the idea. Its outcome (matter must be wave-like) constrains every later step.
2. **39A next.** It is the real mechanism test for H-A, and it tells us whether "matter sets the lapse" holds on our lattice at all.
3. **39C** only if 39A fails or is ambiguous. It is the observational argument against a free N.
4. **39B2** if you want to pursue "particles are lapse patterns" further.

| Outcome | Reading |
|---|---|
| 39A passes, 39C shows scalar radiation | H-A: geometry is fixed by entanglement equilibrium; N should not be a free wave |
| 39A fails, 39B1 passes | Entanglement does not fix geometry here; N as the primary field remains viable, but 39C becomes the decisive test |
| 39B2 finds a discrete mass spectrum | Strong support for H-B: the lapse selects matter |
| 39B2 finds only cutoff-sized boson-star lumps | H-B can hold particles together, but it does not explain particle masses |

## Known limits

- **Scalar sector only.** Everything is in the scalar (Newtonian) sector. Tensor gravity is outside it.
- **Lattice equilibrium is not Lorentz-symmetric.** 39A works on a lattice whose cutoff is not Lorentz symmetric (37A). A mismatch between G_eq and G_ind would be expected, and would not by itself refute A.
- **38B's open Tolman shortfall** (about 0.63 of the prediction) may reappear in 39A as a block-size-dependent residual. The 39A plan should register that possibility.
