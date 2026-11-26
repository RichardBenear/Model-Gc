# EX-GC39C — does a free lapse field radiate? (scalar radiation from an orbiting pair)

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC39C_plan.md) · [design note](../../spec/EX-GC39_design.md). PI CLOSED 2026-11-26. Cloud container, Python 3.13, numpy/scipy. About 0.15 core-h. This is our construction, not LaViolette's Model G.

**Yes. A freely propagating lapse field radiates 1/6 as much power as GR's gravitational waves, so N cannot be a free field with full gravitational coupling.**
- An orbiting pair on the lattice radiates 0.1606 × the GR quadrupole power. The analytic prediction is 1/6, and with the sources' finite size 0.1625.
- If N were all of gravity, binary orbits would decay at 1/6 of the observed rate.
- If N existed alongside GR's tensor waves, the decay would be 16% faster than observed. Binary pulsars agree with GR to about 0.2%, so that is excluded by about 80×.

So the lapse must not radiate on its own. It has to be a constrained quantity fixed by matter (H-A, as in GR), or something must stop its waves from coupling at this strength.

## Results (160³ lattice, h = 1, sponge 20 cells; Gaussians σ = 1, orbit separation 5, Ω = 0.08)

| Case | Lattice power at R = 30 / 40 / 50 | Exact for this source | Lattice ÷ exact | Lattice ÷ GR |
|---|---|---:|---:|---:|
| equal masses (1, 1) | 4.214e-5 / 4.210e-5 / 4.217e-5 | 4.195e-5 | 1.004 | 0.1606 |
| unequal (2, 1), COM fixed | 7.476e-5 / 7.468e-5 / 7.481e-5 | 7.744e-5 | 0.965 | 0.1603 |

- The scalar quadrupole formula, (16/15)Gμ²d⁴Ω⁶, gives 4.369e-5 (equal) and 7.767e-5 (unequal).
- The exact power for the Gaussian source is 0.960 of the quadrupole value (finite size and multipole corrections).
- GR's (32/5)Gμ²d⁴Ω⁶ is 2.621e-4 and 4.660e-4.

**Scoring:**
- **P1 (lattice = exact within 5%): PASS.** 1.004 for the equal-mass pair.
- **P2 (no dipole: unequal case matches the quadrupole formula within 5%): PASS.** 0.962. There is no dipole radiation, as expected for a universal coupling with the COM fixed.
  - The unequal case is 3.5% below its exact value. The exact value includes a small odd-harmonic part; the lattice readout fits only harmonics 2 and 4.
- **P3 (P_s/P_GR = 1/6 within 5%): PASS.** 0.1606 against 0.1667 (3.6% low, of which 2.5% is the finite-size factor).

## Checks and deviations

- **G2: PASS.** The power is the same through the spheres at R = 30, 40 and 50 to within 0.2%.
- **G1: FAIL as run.** The time-domain static field had not settled by t = 200: it was 30–70% below the lattice Poisson solution at r = 25–40. With the sponge, the monopole relaxes slowly.
  - The same lattice operator solved directly (DST Poisson) gives r²·dn/dr = 1.017 GM and 1.005 GM at r = 10 and 20.
  - P1 independently confirms the coupling normalisation.
- **Readout deviations, recorded:**
  - **First run:** the full source was switched on abruptly. The monopole transient swamped the flux; the fit residual was 5× the signal.
  - **Second run:** I removed the orbit-averaged (static, non-radiating) density and ramped the source on over t = 50. That run estimated power from ⟨ṅ²⟩ alone, which is exact only for a pure outgoing wave, and so read 7–18% high at these radii (the near-zone terms 3/(kR)² …).
  - **Final run:** it uses the true time-averaged flux −(1/4πG)⟨ṅ ∂_r n⟩ from fitted harmonics at R ± 1. The ⟨ṅ²⟩ values are kept in the JSON for the record.
  - The physics prediction (1/6) was fixed before any run, and none of these changes touched it.

## What it means (inference)

- **This is the observational argument for H-A over a free H-B.** If the lapse had its own wave dynamics (37C's wave field) with full gravitational coupling, binary pulsars would see it. The lapse must be slaved to matter, as GR's lapse is (it obeys a constraint, not a wave equation).
  - 37A-S2 found the same thing from the other side: the vacuum induces no stiffness for N.
- **37C still holds,** because the moving-source test needs N to *respond* relativistically (retarded, no preferred frame). This result says the responding field must not carry independent energy away. In GR both hold: the radiation is entirely in the tensor modes, and the scalar part is fixed instantaneously in a given slicing.
- **So a substrate in which "the lapse shapes matter" (39B) needs N to be set by the matter and the tensor sector,** not to be a free scalar.

## Files

[Code](../../src/exgc39c.py) · `equal.json`, `unequal.json`, `static.json`
