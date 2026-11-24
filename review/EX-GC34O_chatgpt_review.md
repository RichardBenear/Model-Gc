# ChatGPT review of EX-GC34O

2026-11-24 | Independent review of the exploratory results

**Verdict.** GC34O passes its registered tests and closes the local momentum gap left by GC34N. The scalar packet’s energy, momentum, and stress come from one discrete action and satisfy both of Y1’s local conservation identities to numerical roundoff on the flat background. Proceed to a weak-field sourced-gravity test after fixing the physical packet normalization across spacings.

## Why the result holds up

At each site, the gradient contribution to the normal stress now uses the product of the gradients on the two adjacent links. GC34N used the average of their squares. GC34O extends the one-dimensional correction to three dimensions and shows why the existing shear terms already give the transverse part of the balance. The lapse and shift couplings are unchanged, which preserves the exact energy identity.

The maximum normalized momentum residual across the packet cases is \(2.7\times10^{-15}\). It also vanishes to roundoff on random three-dimensional fields, where the GC34N stress gives 31–41% residuals in the same harness. Direct numerical variation of the action agrees with the analytic sources to about \(10^{-11}\). The random-field control and separate variation check make this stronger than a convergence test on one smooth packet.

## Packet normalization before comparing curvature

The inherited packet generator appears to weaken the physical source on finer grids. It calls a normalized inverse FFT and multiplies by \(N^3(h/L)^3=1\) for \(N=L/h\). For the same resolved Fourier profile in a fixed box, the field amplitude therefore falls by about eightfold when \(h\) halves. The stress is quadratic in the field, and the recorded absolute integrated stresses drop by about 60-fold per halving. The normalized conservation gates in GC34O remain valid. A comparison of gravitational response across grids needs fixed physical source normalization.

## Recommended next experiment

Use the GC34O scalar as a source for Y1’s linear gravity equations. Match the physical packet energy across resolutions and verify the packet width and central wavelength.

1. **Source compatibility.** Check the full forced equations and constraint residuals as the packet moves.
2. **Curvature response.** Compare two preparations with the same initial field profile and different initial field velocities. Calculate their energy, momentum, and stress independently. Compare the difference in gauge-invariant tidal curvature with linearized GR driven by those sources, using matched boundary and initial conditions.

GC34O’s conservation derivation uses the flat-space scalar equations. Keep gravitational feedback on the scalar outside the first sourced test. A self-consistent feedback test, a self-gravitating knot, and the origin of Y1’s geometry remain separate questions.

## Scope and records

I reviewed the plan, action and source code, and committed results. I checked the packet normalization algebraically against the stress records. I did not rerun GC34O.

- [Registered plan](../spec/EX-GC34O_plan.md)
- [Results and verdicts](../results/EX-GC34O/EX-GC34O_note.md)
- [Committed test data](../results/EX-GC34O/results.json)
- [Matter action and source code](../src/exgc34o.py)
- [Inherited packet generator](../src/exgc34m.py)
