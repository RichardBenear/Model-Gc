# EX-GC34Q — kinetic energy and the tidal response on Y1 (EXPLORATORY)

Authorized by Richard in the handoff request. Follow-up to [Claude's handoff](../docs/EX-GC34Q_handoff.md) and the [GC34P review](../review/EX-GC34P_chatgpt_review.md). ChatGPT leads. Cap: 2 core-hours. Report the registered verdict before any interpretation. PI closure is separate.

## Experiment

Use GC34P's normalized massive scalar profile in the same periodic box (L=32, m=1, central k=0.6, width=3). Use both existing profile directions. The standing preparations have identical initial phi and initial velocity c phi, with c=0 and c=1 in the code's units. Each Fourier mode evolves analytically with the flat Klein–Gordon dispersion using its own grid's derivative. Hold the profile and c fixed across spacings. Use h=0.5 and 0.25 and the continuum reference on N=96. No scalar feedback.

Measure the difference, c=1 minus c=0. Because the gravity equations, initializer, and readout are linear, evolve the difference of the two fully evaluated sources. Retain density, momentum and every stress component. Use GC34P's synchronous-gauge minimum-norm constraint initialization, RK4, and gauge-invariant tidal tensor, at t=0, 2, and 4. Both references omit the periodic zero mode.

Use a cube of Fourier modes with integer components at most 14 in magnitude. Audit the omitted squared norm of each nonzero source-difference component separately (excluding zero mode), at t=0, 1, 2, 3, 4, for both directions and all grids. Use the first cut in [14, 18, 22, 25] giving a fraction <=1e-8 everywhere. If none passes, record Q0 fail and stop before tidal scoring.

## Registered gates

- Q0, initial data and coverage. Initial profiles identical within each pair. Total initial momentum in either preparation <=1e-12 of the energy difference, in code units. Independently sum the flat Hamiltonian; the initial energy increment agrees with 0.5 times the integral of phi squared to 1e-12 relative within each grid and with its continuum counterpart to 0.1%. Both flat energies individually drift by <=1e-12 relative over the audit times. Omitted source fraction meets the rule above.
- Q1, compatibility. Initial constraint residual <=1e-10 and maximum <=1e-6 at dt=0.05, both directions and both spacings. For off-axis h=0.5, check dt=0.025 and require an improvement of at least 8x or both residuals below 1e-10. At t=4 the tidal change on halving dt must be below 10% of the spatial discrepancy against the continuum.
- Q2, curvature. Compare the weighted relative L2 error of the tidal difference against continuum linearized GR, both directions at t=0, 2, 4. At h=0.25 each error must be <=2%, and each h=0.5 to 0.25 ratio must be >=3.
- Q3, separate computation. At off-axis h=0.5, use the three modes with greatest final tidal norm. Construct the kernel mode by mode from GC34N's tensor formula, use direct Fourier sums, initialize by a separate least-squares solve and evolve with DOP853 (rtol=1e-10). Require agreement with the main dt=0.05 run to <=1e-3 relative per mode at t=4. Shared dependencies are the scalar profile and 34O source functions and the Y1 tensor equations. This is a numerical cross-check, not an independent replication.

Overall PASS requires Q0–Q3. Report the signed energy increment, curvature-difference norm and numerical scope. This is a check of the response to an assigned stress-energy coupling, not a derivation of gravity or a test of feedback. Register before any computation.
