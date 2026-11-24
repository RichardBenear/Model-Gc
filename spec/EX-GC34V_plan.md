# EX-GC34V — next-order Y1 source compatibility (EXPLORATORY)

2026-11-24 | ChatGPT | Richard authorized the remaining lattice task after GC34U. Cap 2 core-hours.

Test a specific quadratic completion of the linear Y1 equations. The candidate uses the mixed-index continuum Einstein expansion with each derivative replaced by the phase-centered Y1 symbol of the harmonic being differentiated. For the matter source it uses the exact flat GC34O stress and its scalar variation, together with the continuum inverse-metric response of the mixed stress tensor using the Y1 gradient. The scalar response is driven by the GC34R matter variation. No conservation projection is permitted.

Use the single positive-frequency harmonic of a massive scalar with mass 1, amplitude 1 in coefficient units and L=8. Directions are integer modes (1,0,0) and (1,1,1). At orders 0, 1 and 2 the relevant field coefficients are at harmonics 1, 2/3 and 4: the leading metric is at harmonic 2, the scalar correction at harmonic 3, and the next gravitational source at harmonic 4. The calculation is the analytic complex extension of a polynomial interaction coefficient of the real-field equations. In that channel no homogeneous component appears. It does not provide a prescription for the mean density of a real periodic packet.

For h=.5,.25,.125, solve the first metric harmonic in synchronous gauge using Y1, with all ten rows checked. Compute the scalar response at harmonic 3 from the fixed-metric GC34R variation. At harmonic 4 assemble the effective source T1-G2 in mixed indices, with κ=1, and test its discrete divergence using the harmonic-4 derivative. Use a continuum control that replaces K by k and uses the continuum stress and scalar variation. Normalize the divergence by the harmonic-4 Euclidean derivative norm times the sum of the two input tensor norms. Also check the effective source's symmetry after lowering with the flat metric.

Gates:
- V0 first-order solves satisfy all rows within 1e-10, and the continuum control's normalized divergence and symmetry residuals <=1e-10.
- V1 the candidate is compatible in this channel if both normalized residuals <=1e-10 at every tested grid and direction. If V1 fails, report the sizes and refinement ratios, and do not evolve the second-order metric.
- V2 report the Fourier product-rule defect K(2q)-2K(q), and test the geometric quadratic Bianchi combination separately, using the same prescribed first-order metric and the expanded connection products. This distinguishes geometric from matter contributions without making a uniqueness or no-go claim for other completions.

A compatible harmonic channel is a necessary check, not sufficient for real-packet reciprocal gravity. GC34U remains the isolated continuum benchmark.
