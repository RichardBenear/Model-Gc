# EX-GC53 — a link shift for moving matter: can clocks read it, does matter source it, and does it remove 52's preferred frame?

2026-11-27 | EXPLORATORY | Claude | **Plan only. Nothing is run until the PI gives a separate go-ahead.** PI request 2026-11-27 20:48 UTC ("I like it. Write 53."), following the PI-relayed ChatGPT suggestion of a signed link timing offset and a two-ring clock experiment. Cap requested: 2 core-h (expected about 1.3). Our construction; not LaViolette's Model G and not Chain Theory.

## Why

52 found that the 51 model gets the clock part (g₀₀) of a moving mass exactly right. But two bodies moving side by side attract too strongly: R(0.2) = 1.11–1.13 against 1 for Lorentz invariance. The missing piece is GR's shift, or gravitomagnetic part (g₀ᵢ). This is not torsion: GR drags frames with zero torsion. Torsion proper couples to intrinsic spin and is left for a separate study.

53 adds the smallest object that can carry the shift. Each link gets a signed timing offset that is zero unless matter moves. The study then tests it in four steps:
- 53A: can the clocks read it?
- 53B: how does matter couple to it?
- 53C and 53D: with GR's law for it, does a rotating ring give the right clock signal?
- 53E: does the same law remove 52's preferred frame, and with what strength?

**A derivation I tried and could not close.** Before writing this, I tried to derive the shift's strength from the time order in which a moving firing pattern triggers coincidences. It does not close. That ordering is already contained in the time-dependent scalar N, which 51/52 include, and it does not produce a static directional offset. So 53 registers only GR's law (a calibration, not a mechanism). It also measures the strength λ* that Lorentz invariance actually needs in this lattice model. If λ* = 4, the 51 scalar sector is exactly GR-compatible and only the vector part is missing. If not, the scalar sector also needs repair. Deriving λ from event counting stays open.

## The extension (everything else is the 51 model, unchanged)

The extension has three parts: a shift field on links, its coupling to matter, and its coupling to light.

**1. Shift field.** A real variable β_{x,a} on each link (x, x + a), for a = x, y, z. It is odd under link reversal (β_{x+a,−a} = −β_{x,a}), so a change of clock settings at the sites cancels around any closed loop.

**2. Matter coupling** (the ADM form H = N√(m² + p²/…) − β·p, here exact for the 51 metric):
- **The added term.** H_β = −Σ_{x,a} β_{x,a} P_{x,a}, with P_{x,a} = Im(ψ_x†ψ_{x+a}) summed over all four spinor components. This is the lattice momentum carried across the link. The shift acts as a translation, so it is the identity in spinor space.
- **The source.** The matter current is j_{x,a} = −∂⟨H⟩/∂β_{x,a} = P_{x,a}.

**3. Field law (i), the GR benchmark** (harmonic gauge, linear in G). Each link component obeys the same coincidence transport as N:

(1/N̄)∂_t[(1/N̄)∂_t β] = Σ N̄N̄′(β′ − β) − λκ j,

with λ = 4 for GR and κ = 4πG from 51. N̄ is the mean N of the link, and the sum runs over the six parallel neighbouring links. The sign convention is fixed by check C0: the steady β of a uniformly moving point source must equal the boosted Schwarzschild h₀ₓ = −4Gγ²mV/√(γ²X² + ρ²) at linear order.
- Far from a rotating source this gives the gravitomagnetic dipole β = −(λG/2)(J × r)/r³. With λ = 4 that is GR's −2G(J × r)/r³.

**4. Light.** Light obeys the 37B rule with the convective time derivative D_t = ∂_t + β·∇ (central differences along the link direction):

(1/N)D_t[(1/N)D_t u] = Σ N_xN_y(u_y − u_x).
- Eikonal: signals run at N² ∓ β along and against a link, so the crossing times differ by 2β/N⁴.
- The round-trip difference around a closed loop is Δt = 2∮β·dl/N⁴, which ChatGPT's notation writes as 2Σb.

## Parts

**53A. Can a clock read an imposed shift? (1D loop, cheap)**
- **Loop.** A closed chain of 4 × 64 sites (a square path), with N ≡ 1.
- **Pulses.** One site launches a clockwise and a counterclockwise pulse (Gaussian, k₀ = 0.5, width 16) at the same instant and records both returns.
- **Imposed patterns.**
  - (a) Uniform tangential β = 10⁻³.
  - (b) A pure gauge pattern, β = gradient of a random site function χ (amplitude 10⁻³). It must give Δt = 0.
  - (c) β = 10⁻³ on one side of the square only.
- **Weak-field N.** Repeat (a) with N(s) = 1 − 0.01 cos(2πs/256) along the loop.

**53B. Matter's coupling (k-space and real space, cheap).**
- **Packet.** A 52-type packet (m = 0.6, σ = 5, r = 0) at V = 0.2, L = 64.
- **Total current.** Σ P_{x,a} equals ⟨sin k_a⟩ from the FFT.
- **Coupling by finite difference.** ∂⟨H⟩/∂β_{x,a} = −P_{x,a} on 20 random links.
- **Spin current at rest.** A spin-polarized packet at rest carries a circulating current (Gordon decomposition) but zero net momentum. Its P field is reported as the source of the spin–spin shift, a lead into a later torsion study.

**53C. A rotating, spin-balanced ring source (static).**
- **Ring.** A positive-energy Dirac ring with m = 0.6 and r = 0: a Gaussian torus of radius a = 10 and tube width 3, with orbital phase e^{iℓφ} for ℓ = +2, 0 and −2, projected onto the positive band at N = 1.
- **Spin balance.** Spin is balanced by averaging q and P over spin up and spin down along the ring axis.
- **Solve.** n = 129. Solve the static N (51 solver, frozen q) and the static β (law (i) with λ = 4, three DST Poisson solves, Dirichlet faces from the dipole far field).
- **J.** J = Σ r × P, from the lattice.

**53D. Two-ring clock experiment (wave on 1D loop chains threaded through the 53C fields).**
- **Enclosing loop.** A stationary square loop concentric with the ring in its plane, with half-widths 20, 28 and 36. Clockwise and counterclockwise pulses are timed by one clock.
- **Small probe loop.** A fixed 8 × 8 square:
  - in the ring's plane at distance d = 16, 24, 32 and 40 along x, with its normal along the ring axis;
  - on the axis at heights d = 16, 24, 32 and 40, with its normal along the axis.
- **Rotated apparatus.** The ring tilted so its axis is along (1, 1, 1)/√3, with one enclosing loop of half-width 28 in its plane (loop sites rounded to the lattice).

**53E. 52D with the shift (full dynamics).**
- **Setup.** The 52D setup unchanged: n = 65, m = 0.6, σ = 5, r = 0, transverse separation 14, T = 16, dt = 0.1.
- **Fields.** N and the three β components are all dynamic, starting from the steady comoving solution of the frozen rigid pair (N as in 52, β from law (i) with ∂_t → −V∂_x).
- **Runs.** V = 0, 0.1, 0.2 and 0.3, with λ = 4. Each run is a pair run plus a single-packet run for the self-force subtraction, as in 52.
- **Spin.** Both packets have the same spin, so the 52 side shift leaves their separation unchanged. At V = 0.2 a spin-down repeat gives the spin average.

## Predictions (scored)

- **P1 (readout, 53A).**
  - Pattern (a): Δt_wave = 2Σβ/N⁴ within 2%, for N ≡ 1 and for the weak-field N.
  - Pattern (b): |Δt| < 10⁻³ of pattern (a)'s value.
  - Pattern (c): Δt equals one quarter of (a)'s value within 2%.
- **P2 (coupling, 53B).**
  - Σ P = ⟨sin k⟩ to 1e-10.
  - The finite-difference ∂⟨H⟩/∂β equals −P to 1e-8 on every tested link.
- **P3 (rotating ring, 53C/D, λ = 4).**
  - (a) Δt(ℓ = 0) = 0 to 1e-12 (absolute, in units of Δt(ℓ = 2)), and Δt(+2) = −Δt(−2) to 1e-6 relative. The scalar N must be identical for ±ℓ to 1e-10.
  - (b) Small loop: the log-slope of |Δt| against d over 16–40 is −3 ± 0.15, in the plane and on the axis.
  - (c) The ratio of the axis signal to the in-plane signal at equal d is −2 ± 0.1 (the dipole pattern).
  - (d) At d ≥ 24, and for the enclosing loops, Δt equals the analytic dipole value 2∮β_dip·dl (β_dip = −2G(J × r)/r³, with the lattice J) within 5%.
  - (e) Tilted ring: Δt equals the aligned value scaled by the loop's actual enclosed geometry within 3%.
- **P4 (preferred frame removed, 53E).**
  - R_dyn(V) = γF₁(V)/F₁(0) = 1 within 2% at V = 0.1, 0.2 and 0.3.
  - At V = 0.2 the spin-up and spin-down R_dyn agree within 1%.
- **P5 (the strength Lorentz invariance needs).** R depends linearly on λ. Combine 52's R_dyn(0.2) (λ = 0) with 53E's (λ = 4) to get λ* where R = 1. The prediction is λ* = 4 ± 0.3.

## Checks (separate paths)

- **C0 (sign and normalization, before the scored runs).** The 53E steady β of a single point-like moving source at V = 0.2, n = 65, equals the boosted harmonic h₀ₓ above within 3% at transverse distances 8–24.
- **C1.** The 53C β solve is checked against an independent conjugate-gradient Poisson solve of one component to 1e-8 relative. I learned in 52 to keep κ large enough that the floating-point floor sits well below this.
- **C2.** 53E at dt/2 for V = 0.2: F₁ changes by less than 1%, with the same κ (52's C3 bug is fixed).
- **C3.** With β ≡ 0, 53E's code reproduces 52's F₁ at V = 0.2 to 1e-10.

## What the outcomes mean

- **P1–P4 pass and λ* ≈ 4.**
  - The 51 scalar sector plus a standard shift gives Lorentz-invariant gravity for moving matter to first post-Newtonian order. The clocks read frame dragging as a one-clock Sagnac delay, the delay flips with the source's rotation, and it falls off as a dipole.
  - The open question is then one only: what in the event picture produces the shift with λ = 4. That is a derivation problem for the next step, possibly in a non-symmetric rewriting rule.
- **λ* ≠ 4.** The scalar sector itself is not GR-compatible for moving sources (for example through its transverse spatial part). I would locate the mismatch before adding anything else.
- **P3 fails while P4 passes, or the reverse.** The coupling or the readout is inconsistent between translation and rotation. I would trace that first.

## Outputs

`src/exgc53.py`, `src/exgc53_recompute.py`, `results/EX-GC53/`, note, LOG row and HANDOFF section.

## Amendment A1 (2026-11-27, after 53A/53B/C0 had run and before 53C–E, C1–C3)

**What had already run.** 53A, 53B and C0 were run first, with the final code. 53A gave registered fails:
- pattern (a): 0.983 of the readout (pass);
- pattern (b): 0.76% of (a) (fail; the limit was 0.1%);
- pattern (c): 0.904 of the readout (fail);
- weak-field N: 0.825 of the readout (fail).

**What diagnosis showed** (recorded as post hoc in the note):
- *The weak-N miss* is an artefact of the clock's centroid density. u² + w² is not the conserved energy when N varies. With the conserved energy density, the weak-N ratio equals the uniform one (0.938 against 0.936).
- *The remaining 6%* is lattice dispersion at k₀ = 0.5. The lattice eikonal gives a factor cos k₀ / cos²(k₀/2) = 0.935, and the wave matches it to 0.1%.
- *So the clock does read Σβ*, through a known lattice factor. The registered numbers stand.

**Changes for the parts not yet run:**
1. **53D small loops.** These are scored with the readout sum 2Σβ/N⁴ (the quantity 53A validated up to the lattice factor), not with a wave. An 8 × 8 loop has a perimeter of 32 sites, too short for 16-site-wide pulses.
2. **53D enclosing loops.** P3(d) is scored on the readout sum against the analytic dipole. The wave delay, now with the conserved-energy centroid, is reported against readout × 0.935 but not scored.
3. **53D tilted ring, P3(e).** The ratio (readout sum / analytic dipole sum along the same staircase loop) for the tilted ring equals that for the aligned ring, along a staircase loop of the same radius, within 3%.
4. **C2.** Shortened to fit the cap. 53E at V = 0.2 runs to T = 8 at dt and at dt/2, with the same κ. The packet-1 momentum difference P₁(t) at shared times must agree within 1% at t ≥ 2.
5. **C3.** With the shift off, the P₁(t) values of 53E's V = 0.2 setup over 0 ≤ t ≤ 2 equal 52's recorded P₁(t) to 1e-10.
