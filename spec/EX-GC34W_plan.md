# EX-GC34W — cubic lattice action and vacuum conservation (EXPLORATORY)

2026-11-24 | ChatGPT | PI: “go ahead with the experiment”. Cap 2 core-hours.

Keep the GC34N linear Y1 operator and the GC34O scalar coupling. Start with the necessary pure-vacuum part of the proposed action test. A failure with scalar field identically zero rules out this candidate combined completion before matter or time evolution is tested.

Define one short-range cubic gravity candidate from the Einstein–Hilbert connection-product action after dropping its boundary term. Use Y1's staggered one-cell derivative for each spatial derivative, and average each field or derivative to the vertex of its cell before multiplying the three factors. The placements follow index parity; an average in direction a contributes cos(q_a/2) on a position-referred Fourier amplitude. In constructing the cubic coefficient, expand the inverse metric, volume factor and connection through the required order. Reconstruct the source by varying the ten independent metric amplitudes. Use equal lattice and continuum conventions (overall action normalization is immaterial to the test).

Put two linear vacuum waves on shell (each has omega=|K| and transverse-traceless spatial polarization). Obtain the quadratic source at the sum wavevector and frequency from the action's three-leg coefficient. Test its contraction with all four linear gauge directions. On the two free wave solutions, the first nonlinear gauge variation of the quadratic action vanishes, so a compatible cubic action must give zero. One channel repeats the axis or diagonal wave to test the self-harmonic limit. A second uses nonparallel directions (1,0,0) and (0,1,1); use both polarizations of each wave. Set L=8 and h=.5,.25,.125. Compare to the continuum control, with no averaging and K=k.

- W0: both input waves satisfy the full linear kernel to <=1e-12 relative; continuum on-shell cubic contractions <=1e-10 normalized.
- W1: action-derived source is symmetric by construction; require finite-spacing gauge contractions <=1e-10 in every channel. Normalize by the source norm times each gauge vector norm, and report the source norm.
- W2: if W1 passes, proceed to the mixed matter test before time evolution. If W1 fails, report the obstructing channels and grid scaling. Do not project the source or proceed to a reciprocal run. Use a direct evaluation with a gauge test leg as a separate check of the variational contraction.

This tests one specified nearest-neighbor interpolation of the cubic action, not all possible actions.
