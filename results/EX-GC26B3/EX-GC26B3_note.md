# EX-GC26B3 — ramped core-force calibration

Completed 2026-11-20 (PI local date). Four cases, 1.331 measured core-hours (<2 cap), unchanged B2 spatial stencil and hash-verified rest states. Plan registered locally before execution; subsequently checkpointed to GitHub. No precision gate. This is the massless Skyrme reference model, not Model Gc.

## Primary independent measurements

The squared-deficit force was smoothly ramped over t=0–1 and held thereafter, with F(t)=f(t)Q(t), Q=integral (1-phi0)^2. Initial nominal F/E0=.025. RK4 evaluates force at each stage. Paired no-force controls at both resolutions were subtracted; raw and corrected results are effectively identical. Impulse is trapezoidal integration of saved samples. Late velocity is a linear fit to the global energy centroid for t>=4; momentum averaged over those identical samples. Units c=1.

| Quantity | h=.2 | h=.15 |
|---|---:|---:|
| E0 | 1.20240177 | 1.21710441 |
| T | 5.0000 | 5.0025 |
| Applied impulse/E0 | .11498644 | .11494278 |
| Gained P/E0 | .11094992 | .11282235 |
| Gained P / applied impulse | .96489570 | .98155228 |
| Impulse shortfall | 3.5104% | 1.8448% |
| Late speed | .09749687 | .10145371 |
| P/v / E0 | 1.00605807 | 1.00266086 |
| P/(gamma v) / E0 | 1.00126504 | .99748739 |
| Acceleration fit t>=1 | .02452622 | .02511639 |
| Force/acceleration / E0 (t>=1) | 1.04229405 | 1.01801730 |
| Acceleration fit t>=2 | .02493996 | .02542395 |
| Force/acceleration / E0 (t>=2) | 1.03366161 | 1.01379757 |

## Error budget for comparative 26E

- Impulse transfer differs by 1.67 percentage points across grids and improves with refinement; the shortfall ratio is .526 vs the h^2 ratio .5625. Compatible with leading discretization error, not a convergence proof.
- Direct mass differs by 0.34 percentage points between grids (gamma-corrected, 0.38 pp). This calibration is for a single B=−1 state and does not guarantee the same error for other species.
- Force/acceleration mass is 4.23%/1.80% above E0, whereas direct P/v mass is only .61%/.27% above. Do not conflate the two. Fit-window sensitivity: acceleration increases 1.69%/1.22% (coarse/fine) when start moved from 1 to 2.
- Regional late speeds (R≈3,5) relative to full centroid: coarse −3.17%,−1.10%; fine −3.23%,−1.22%. Use the same global centroid convention for comparisons; treat region sensitivity as a diagnostic, not an uncertainty proven to be common across species.
- Bare energy rises .8408%/.8524%; the source-inclusive ledger H+V+integral f'(t)D(t)dt closes within 1.07e-6/8.46e-7 E0 (sampled quadrature). Null energy drift <4e-13; null momentum/E0 <1.5e-8, transverse forced momentum/E0 <1.5e-8. All initial/final geometric degrees −1 in all three preimage checks; zero bad tetrahedra.

Recommendation: proceed with comparative 26E, with direct mass and momentum-transfer readouts separate, same normalization and fit windows, and a coarse/fine check on representative states. Treat a few-percent force-response systematic as the observed baseline, not as a universal bound or a new precision gate. No further stencil refinement at this stage. Any apparent species-dependent effect comparable to these errors should be treated as unresolved. Stronger differences could justify targeted follow-up. This does not establish Model Gc inertia or equivalence; it calibrates the reference solver.

Reproduce `python src/exgc26b3_summary.py` and `python src/exgc26b3_figure.py` from the four JSON records. See summary.json for full numerical metrics, per-run hashes and raw sample histories in individual records. Figure: paper/figures/exgc26b3_calibration.svg.
