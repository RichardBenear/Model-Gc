# EX-GC35M — local constraint damping on the Y1 stencil

2026-11-24 | EXPLORATORY | ChatGPT | [Registered plan](../../spec/EX-GC35M_plan.md). **PI CLOSED 2026-11-24.**

**Scalar residual and tidal gates failed. The fixed full-stress knot sources also failed.** Damping reduced the free scalar defect from GC35C but left sustained late growth over the tested window. The approximate stress coupling was retained. No divergence projection or moving-knot run was used.

## Scalar controls

The synchronous metric RHS adds +α trace-basis times the Hamiltonian residual. The auxiliary velocity RHS adds −β times the tensor-weighted adjoint of the momentum constraint operator applied to the residual. Both α and β are .12. In Fourier notation this changes ḣ to v+α I r₀ and v̇ to a−β W⁻¹RB†r. The two correction operators are a local trace update and a first-difference adjoint. There is no per-mode inverse or reset after initial data. The scalar benchmark still uses the same band-limited evaluated source and RK4 as the prior experiments. The trace and adjoint updates have TT content below 1.34e−15 in the direct mode audit. The constraint acceleration block is zero within 3.56e−15. A Brillouin-cube sample gave maximum RK4 amplification 1.000000 and 1.000000 at Δt=.05. This is an eigenvalue screen rather than an all-mode nonlinear stability result.

The reported constraint uses the actual corrected ḣ rather than the auxiliary v. The tidal tensor uses second differences of the metric positions, including the damping update.

| h | Peak original residual | Residual at t=20 | Reconstructed tidal error at t=4 | At t=20 |
|---:|---:|---:|---:|---:|
| .5 | 0.1732 | 0.1732 | 13.5% | 108.8% |
| .25 | 0.04494 | 0.04494 | 3.38% | 27.6% |

Peak residual refinement was 3.85×, while the finer peak remained above 10⁻³. The h=.25 reconstructed tidal error at t=20 improved from 398% in the matched free replay to 27.6%. The h=.125 projected reference was 2.25% at this time, and even the matched-spacing projected reference was 8.94%. Thus the damping accuracy gate fails independently of the spacing mismatch with GC35G.

## Held-static knot sources

The h=.2 GC35H candidates were regenerated to the same tight force tolerance. Their full stress and energy-only controls were frozen, evolved in synchronous variables, and compared using the original constraints through t=20. This is a compatibility audit with fixed densities, not an evolution of either knot. These residuals divide by the original source norm without the small scalar experiment's perturbative scaling.

| Knot | Full source at t=5 | At t=20 | Peak energy-only control |
|---|---:|---:|---:|
| B1 | .0915 | .1843 | 3.17e−14 |
| B2 | .0662 | .1489 | 2.89e−14 |

The static audit shows that this damping choice cannot make the original equations compatible with the full held stress. The passing energy-only control qualifies the diagnostic, and the failed scalar and full-source gates leave moving knots gated. The run tests a specific local correction, with no conclusion about other local formulations.

Recorded CPU use was 553 seconds. [Source](../../src/exgc35m.py) · [Static driver](../../src/exgc35m_static.py) · [Coarse](n24.json) · [Fine](n48.json) · [Static B1](static_b1.json) · [Static B2](static_b2.json)
