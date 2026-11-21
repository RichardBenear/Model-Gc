# EX-GC26B3 — comparative force calibration

PI approved 2026-11-20. No precision gate; characterize systematic errors before comparative EX-GC26E.

Reuse hash-verified B2 mirror-symmetric relaxed states and dynamics, c=1, L=7, h=0.2 and 0.15. Four cases: force and matched no-force at each resolution, initially at rest, T=5 (nearest step), dt=0.05h. Only squared-deficit coupling V=-f(t) integral x(1-phi0)^2, f(t)=0.025 E0/Q0 times half-cosine ramp from zero to one over t=1. Evaluate f at all RK4 stage times. Sample approximately every 0.1. Hard aggregate budget 7200 CPU seconds; per-case limits 1050s coarse and 2450s fine (7000s total), leaving 200s overhead. No reruns without budget check.

Primary independent readouts: endpoint delta P_x / integral F dt, F=f(t)Q(t); and P_x/v_x (also gamma-corrected), with v from linear fit to energy centroid over t=4–5 and P averaged over the same window. Subtract matched control drift for both; report raw values too. Fit quadratic centroid and linear momentum only post-ramp (t>=1), compare t>=2 and regional centroid fits. Monitor degree, transverse momentum, regional centroid dependence, bare H and source-inclusive energy ledger H+V plus integral f'(t) integral x(1-phi0)^2 dt. Report numerical error budget and implications for relative measurements; no claim of new physics or precision qualification.
