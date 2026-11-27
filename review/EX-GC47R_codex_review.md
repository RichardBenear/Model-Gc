# Independent review of EX-GC47R

**GC47R confirms inward motion after correcting the release packet and the microscopic strain coupling. The six-site force comparison needs a narrower accuracy claim and a sign correction in the saved continuum diagnostic.** [Plan](../spec/EX-GC47R_plan.md) · [Results](../results/EX-GC47R/EX-GC47R_note.md)

## What survives the corrections

Building and tracking the packet about its own center removes the artificial jump identified in GC47. Across the tested widths, separations, windows and box sizes, inward acceleration agrees with the effective pair-band prediction to within 4.4%. B1 accelerates at 90.2–94.6% of B2's acceleration under matched conditions. The corrected runs support a composition difference for these two model bodies. They share the same effective Hamiltonian construction and are sensitivity checks, not independent physical measurements. The acceleration is about 1.15–1.32 × 10⁻⁴ of the chosen lapse-gradient reference. Interpreting that absolute scale as a failure of relativistic inertia still requires a separately calibrated signal speed. [Release results](../results/EX-GC47R/EX-GC47R_note.md)

Including the derivative of the band projector brings the zero-momentum strain charge into agreement with the uniform energy derivative to better than four parts in ten billion in the reported cases. The change in the L = 32 static interaction is at most 0.017%. This resolves the earlier normalization failure and shows why its effect on the attraction was small. [Form-factor data](../results/EX-GC47R/F.json)

The registered energy-drift gate failed. The worst reported error is 2.3 × 10⁻⁷ of the potential-energy scale, versus a 10⁻¹⁰ gate. The code measures final minus initial energy and compares endpoint displacement when halving the time step. Those are encouraging diagnostics but do not bound the largest energy excursion during a run or directly test convergence of the fitted acceleration. Record both quantities in the next dynamics run. [Release implementation](../src/exgc47r.py)

## What six lattice spacings actually resolve

For B1 against B1 at L = 48 and r = 6, the saved force is inward. Its magnitude differs from that of a point source in the same periodic lattice by **0.0083%**. It differs from the central difference of the continuum potential −C/r by **1.89%**, and from the analytic inverse-square magnitude C/r² by **4.80%**, where C = G_eff(E² + 2EQ_a). All three numbers describe the same saved force. The distinction arises because the unit-step central difference gives C/(r² − 1). At r = 8 the discrepancy with the analytic inverse square falls to **0.84%**. These comparisons quantify numerical behavior in lattice units. The calculation uses a prescribed lattice Poisson kernel. [Force data](../results/EX-GC47R/force.json) · [Force implementation](../src/exgc47r.py)

The code assigns the continuum comparator the opposite sign from the inward force. Consequently the saved signed ratio at six sites is −1.01890, while the note quotes its positive magnitude. Correct the comparator or label the comparison as a magnitude. The microscopic-to-periodic ratio uses consistent signs and its registered pass stands. Increasing box size removes periodic-image effects while leaving the lattice spacing fixed. At r = 6 the reported magnitude ratio to the differenced continuum crosses unity as the box grows. A near-perfect match in one box can therefore reflect cancellation between finite-box and discretization errors. The x-axis sweep alone cannot establish direction-independent continuum accuracy. [Saved forces](../results/EX-GC47R/force.json)

## Next decision

I would close the two GC47 correction items, retain the energy-drift failure on the record, and carry the force-reporting fix into the next analysis. The next physics test should measure the signal speed independently in the same distance and time units as the massive pair dispersion. Compare c_signal² with 2 t_pair E_rest for both B1 and B2 before any adjustment to G. A common mismatch would identify a matter–signal calibration problem. A difference between the two ratios would survive a common change of units. If short-distance gravity becomes the focus, vary lattice spacing at fixed physical separation and compare axial and diagonal forces.

I reviewed the plan, results and implementation and recalculated the force comparisons from the saved output. I did not rerun the simulations.
