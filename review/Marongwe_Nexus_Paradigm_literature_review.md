# Literature review: Marongwe's Nexus Paradigm (NP) and its relevance to Model Gc

Claude, 2026-11-23. Requested by the PI (Richard): "Take a look at papers by Marongwe / Nexus-Paradigm (NP) /
Quantum Spatio-Temporal Dynamics (QSD) framework and see if there are any ideas that may be applicable to our work."
Literature only: no runs, no code changes, no compute. Not a registered study.

## Conclusion

**Nothing in the Nexus Paradigm addresses our current blocker, and I found no idea in it worth a test.** NP has no
lattice dynamics, no lapse or shift, no Hamiltonian constraint and no gravity kinetic term beyond writing
√g (R − 2Λₙ). Its "quantisation" is a discrete label n on Einstein-manifold solutions (Λₙ = n²Λ), not a quantised
geometry with its own dynamics. The parts that touch our past lessons go the wrong way: its graviton mass would spoil
the two-polarisation, γ = 1 results that EX-GC34C–34E were built to keep, and its light bending equals GR only because
the metric is matched to Schwarzschild by hand. The work has little independent uptake (about 14 non-self citations
across 14 INSPIRE records) and I found no independent technical review or follow-up by other groups.

## What was read

| Source | Venue / status |
|---|---|
| A Covariant Canonical Quantization of General Relativity (2018) | Adv. High Energy Phys. 2018, 4537058 (Hindawi), full PDF read |
| Energy-Momentum Localization in Quantum Gravity (2024) | Phys. Scripta 99, 025306; preprint versions arXiv:2402.10910, Preprints.org 202311.1535 (first half verbatim) |
| Turning Galaxy Rotation Curves into Radial Cosmic Chronometers (with S. Kauffman) | arXiv:2604.17597; Semantic Scholar lists venue *Galaxies* (2026) |
| SAGS 2025 talk "Quantum Gravity Outperforms General Relativity" (with Corda, Tshipa) | abstract only |
| Publication list | academia.edu profile, INSPIRE API, Semantic Scholar API |

Not retrieved: the 2021 AHEP paper (Wiley page refused), ResearchGate pages. The 2014/2015 IJGMMP papers ("The Nexus
graviton", "The Schwarzschild solution to the Nexus graviton field") were seen only through later papers that restate
them. Equations from the 2018 PDF come from a poor text extraction; I checked them for consistency with the 2024
preprint, which restates the same set.

"QSD" as a search term mostly returns an unrelated field (quantum spatio-temporal dynamics in tensor networks,
turbulence, open systems). Marongwe's own label is "Nexus Paradigm"; I found no separate QSD body of work by him.

## Peer review and citations

- Journals: IJGMMP (2014, 2015), IJMPD (2017, 2023), AHEP (2018, 2021), Physica Scripta (2024), MDPI *Universe*
  (2025) and *Galaxies* (2026), Int. J. Astron. Astrophys. (SCIRP, 2013). Several recent items are arXiv-only.
- INSPIRE (fetched 2026-11-23): 14 records, 34 citations in total, 14 without self-citations. Highest: the 2015
  Schwarzschild paper (7, 2 non-self) and the 2017 IJMPD paper (5, 4 non-self). Semantic Scholar: 21 papers, 29
  citations, h-index 3.
- I did not check who the 14 non-self citers are. Regular co-authors include S. Kauffman, C. Corda, M. Tshipa and
  M. Mafu. I found no paper by another group that tests, extends or critiques NP.
- So: journal-published, but essentially unreviewed by the field beyond the original referees.

## The construction, briefly

1. Displacement "square root": requiring Δs² = η_μν Δx^μ Δx^ν to have no cross terms gives Δx^μ = a γ^μ with
   {γ^μ, γ^ν} = 2η^μν (Dirac matrices). This is a relabelling, not a derivation of a physical spinor geometry.
2. "Nexus graviton": a sinc-windowed Fourier packet of 4-space with k = nπ/r_H, n = 1 … ~10⁶⁰ (Hubble radius over
   Planck length). Called Bloch states, but no periodic potential is identified. Spin 2 is claimed by adding four
   spin-½ objects, which is not a valid angular-momentum argument.
3. Vacuum as a Ricci soliton / Einstein manifold: Ricci flow ∂ₜg = −2 Ric, rewritten ∂ₜg = −(c r_H/3) G, and the
   field equation G⁽ⁿ⁾ = κT + (n² − 1)Λ g⁽ⁿ⁾. The n-dependence is algebraic, so different "states" are different
   cosmological constants, not different geometries evolving in time.
4. Canonical quantisation is formal: [q̂, p̂] = iħ on the label variables, Hamiltonian density set equal to the
   Lagrangian density "because the graviton does not translate". No Hamiltonian constraint, no lapse/shift, no gauge
   (re-slicing) analysis.
5. The 2024 paper adds a Helmholtz equation ħ²□g + m_G²c² g = 0 with a graviton mass m_G ~ hH₀/c².
6. Phenomenology: matching n² = c²r/GM to Schwarzschild gives orbits r_n = n² GM/c², deflection θ = 4/n² (= 4GM/c²r),
   and a critical acceleration a₀ = cH₀/2π ≈ 1.05e-10 m/s² (MOND-like Tully-Fisher). Later papers use a₀ and an
   "evolving" BTFR to date galaxies, and claim EHT fits favouring NP over GR (Bayes factor ~10³⁶).

## Problems I found (critical reading)

- **Internal contradictions (2018 paper).** It equates the Einstein tensor with the Ricci tensor (needs R = 0) while
  its own field equation G = Λg gives R = −4Λ. The ground state n = 1 has zero effective Λ, yet is called the source
  of dark energy. The horizon is claimed at half the Schwarzschild radius, but the paper's metric factor 1 − 2/n²
  vanishes at r = 2GM/c², the usual radius. ISCO is put at n = 1 (r = GM/c²), inside the claimed horizon and far from
  GR's 6GM/c².
- **Light bending is input, not output.** θ = 4GM/c²r follows from the Schwarzschild matching; no PPN γ is computed.
- **Massive graviton.** A Fierz–Pauli-type mass gives five polarisations and the vDVZ discontinuity (in linear theory,
  light bending stays 3/4 of GR as the mass goes to zero) unless screened; a generic mass term also has the Boulware–Deser ghost
  nonlinearly. None of this is discussed. (Standard results, my background knowledge, not from the papers.)
- **No dynamics of geometry.** Ricci flow is first order in time and parabolic: it smooths metrics, it does not
  propagate waves. It cannot serve as a gravity kinetic term.

## Relevance to Model Gc, item by item

Our blocker (EX-GC34E–34H): a lattice gravity kinetic term that (i) respects re-slicing (no growing re-slicing mode,
which the per-tet DeWitt term has), (ii) gives every metric edge combination dynamics (cube-mean B leaves one
non-kinetic body-diagonal combination with a growing Jordan chain), (iii) keeps two light-speed modes, and (iv) keeps
the 34C γ = 1 conditions (pure clock plus a single covariant gravity energy).

| NP idea | Could it help? | Verdict |
|---|---|---|
| Discrete n, Λₙ = n²Λ | Global label on solutions; no local lattice, no edges, no kinetic term | Not applicable |
| Ricci flow as the "dynamics" of the vacuum | First order, diffusive; would add damping, not a hyperbolic kinetic term. It does not touch the re-slicing symmetry | Not applicable; would break wave propagation |
| Massive graviton m_G | A mass term would give the non-kinetic 7th-edge combination a "frequency", but by breaking the re-slicing gauge symmetry outright: extra polarisations, γ ≠ 1 (vDVZ), contrary to 34C/34E | Rejected; it hides the 34H problem rather than solving it |
| Δx^μ = aγ^μ (Clifford square root of the line element) | Superficially like using tetrads instead of metric edges. A tetrad / first-order (Palatini, Regge–Cartan) lattice formulation is a real option in the lattice-gravity literature, but NP does not construct one and gives nothing to reuse | Not from NP; see "pointer" below |
| a₀ = cH₀/2π, evolving BTFR, galaxy ages | Cosmological phenomenology; Model Gc has no cosmology yet, and these fits are not tests of a lattice rule | Out of scope |
| EHT claims (horizon at ½ r_s) | Contradicted by the paper's own metric | No |

**Pointer, not from NP.** The only adjacent direction the reading prompted is already on our list in a better form:
route 2 of 34H ("derive the kinetic term from an exact lattice re-slicing identity"). If we want a different variable
set than route 1 (drop the body-diagonal edge), the established option is a first-order lattice formulation with
tetrads/frames and a connection, where re-slicing and local Lorentz invariance are built into the variables. That would
need its own literature check (e.g. lattice Palatini / Regge–Cartan / spin-foam-adjacent work) and its own plan; it
owes nothing to NP.

## Recommendation

Do not spend compute on NP-derived ideas. Keep the 34H decision between (1) a six-metric-component variable set and
(2) a kinetic term from an exact lattice re-slicing identity. If useful, a short paper-only screen of first-order
(tetrad plus connection) lattice gravity could be a separate, PI-approved exploration.

Sources: INSPIRE literature API (query "a Marongwe"); Semantic Scholar author 102588570; inspirehep.net file
813044122f3ec91bd3ebb38a7b326419 (AHEP 2018 PDF); alphaxiv.org/abs/2402.10910; preprints.org/manuscript/202311.1535/v2;
arxiv.org/pdf/2604.17597; indico.global/event/13404/contributions/135807; independent.academia.edu/StuartMarongwe.
