# EX-GC34S — constraint audit before reciprocal evolution (EXPLORATORY)

PI authorization: “go ahead with the remaining step”. ChatGPT. Cap: 2 core-hours.

GC34R evolved a scalar tangent on the metric sourced by the flat scalar. The outstanding step is to feed the scalar response back into gravity. The current repository supplies a linear vacuum operator and a matter action through first order in the metric. This registration tests the direct resourcing option before treating that option as reciprocal gravity.

## Method

Use the GC34R standing packet, lattice sizes N=32 and 64, and registered Y1 metric and scalar tangent equations. Run to t=2 with dt=.025, with a dt=.0125 check on N=32. At t=0,1,2 compute the flat action’s source vector and its directional derivative along the scalar tangent. Derive the time derivative by polarization using the analytic flat acceleration and the tangent acceleration. Apply GC34N's local energy and momentum identities. Separately compute the identities on the flat field as a control.

For amplitudes epsilon=1e-4 and 5e-5 evaluate the complete flat-source formula on the corrected field and compare the balance defect with epsilon times the tangent result. Before interpreting a nonzero result, confirm the sign and scale for the source conservation identity independently from the exact discrete energy balance in terms of the scalar equation residual and the site average of neighboring residuals where appropriate.

## Gates and decisions

- S0: flat source identity residuals <=1e-10 at both spacings and all sample times. The metric constraint residual remains <=1e-6.
- S1: report whether the resourced tangent satisfies the original linear identities. A normalized tangent defect above 1e-4 at either positive time on the fine grid, with the timestep check agreeing to 1% on the coarse grid for defects above 1e-4, rules out direct resourcing for this setup. Otherwise the result is unresolved unless every normalized defect is <=1e-10.
- S2: the defect from the corrected scalar at epsilon=5e-5, divided by epsilon, agrees with the tangent defect to <=1% relative for both identities at positive times with a measurable defect. Record both amplitudes and sign. Separately check the energy defect against the analytic local product of field velocity and the Klein–Gordon acceleration residual.
- If direct resourcing is ruled out, do not project away the defect and report that as physical feedback. Record the missing second-order gravity and metric-dependent source terms and the periodic homogeneous-mode issue.

The completed test is a diagnostic of an actual resourcing prescription. It does not establish reciprocal nonlinear evolution if the consistency gate fails.
