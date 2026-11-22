# EX-GC29A — Natural particle response cutoff

Exploratory; 2026-11-22; design/run ChatGPT. PI authorized. Reviewed (review/EX-GC28H_29A_claude_review.md); **PI CLOSED 2026-11-22 (Richard: "go ahead and close 28 and 29A and we will defer 29B")**. **Positive for a size-dependent elastic roll-off; negative for a demonstrated universal frequency cutoff.** No claim of a bulk ZPF regulator or induced Newton constant.

## Method and controls

See spec/EX-GC29A_plan.md (registered before production) and spec/EX-GC29A_dispersion_check.md (registered after initial results, before diagnostic). Twenty cases: C=0 and C=0.17839644285107545, R=12/24, h=.1/.05/.025/.0125, plus two comoving dilations per C at fixed dimensionless shape. The GC28A canonical radial Hessian is used; positive eigenvalues, spectrum reconstruction, uniform-probe Parseval and quadratic-form assertions passed. Elastic radial Fourier–Bessel transform for baryon, energy, P1=1−cos f and P2=P1², with matched RMS Gaussian. Inelastic linear scalar probe is either spatially uniform or j0(omega*r) on-shell; the canonical localized drive is J=r² P'(f)/sqrt(6 pi K). Weights h|v dot J|² (or with j0) and cumulative moments M_p=sum w omega^p, p=0,2,4. Matched l=1 vacuum Hessian with the same J is a geometric control, not an actual vacuum defect. Neither these moments nor their difference are physical ZPF energies; no EM coupling or higher angular channels included.

## Elastic results (R=24, h=.0125)

| Operator | C=0 k50 | C=0 RMS | massive k50 | massive RMS | massive k90 |
|---|---:|---:|---:|---:|---:|
| Baryon | 1.4013 | 1.0575 | 1.4946 | .9837 | 2.8172 |
| P1 | .4972 | 3.9585 | .7626 | 2.1414 | 1.6412 |
| P2 | 1.3703 | 1.0950 | 1.5029 | .9799 | 2.8342 |
| Energy | 1.2041 | 1.4254 | 1.2982 | 1.2359 | 2.6650 |

k50 and k90 denote first half-power and 90%-suppression crossings in dimensionless momentum. At k=8 all massive form-factor powers are below 1.3e-5. Massive P2 k50 under 1.5 enlargement / baseline / 1.5 shrinkage is 1.0019/1.5029/2.2544; k50 times RMS invariant under comoving dilation. This is a scale covariance check, not proof a new independently generated physical length exists. The matched non-topological Gaussian predicts k50=sqrt(3 ln 2)/RMS, similar for P2 and baryon. Long-tailed massless P1 is box sensitive (k50 .5882 at R12 vs .4972 at R24), and its infinite-domain RMS is not finite for the massless f~r^-2 tail. No universal operator-independent scale. The Skyrme energy balance sets a size via model couplings; topology alone does not fix it.

## Radial dynamical results

Massive P2 uniform probe: M0=.0326931, M2=.0624213, M4=.246109 at Omega=32; essentially saturated by 16, with 90%-cumulative frequencies 1.993, 2.649, 5.144. P1 uniform corresponding 90% frequencies .970, 3.043, 3.962. Initial on-shell moments did not saturate in higher powers; coarse lattice phase mismatch can mimic high-frequency coupling despite small relative eigenfrequency error. Additional diagnostic: cubic interpolate fixed h=.0125 profile and smooth defect potential, keep exact centrifugal term; h=.0125/.00625/.003125, eigenmodes below 64, compare incident k=omega with lattice-dispersion matched k=2 asin(h omega/2)/h. This is a fixed-background refinement, not a new BVP.

| Massive P2 M4 | Omega=8 | 16 | 32 | 64 |
|---|---:|---:|---:|---:|
| h=.0125 continuum defect | 1.126 | 3.772 | 14.777 | 431.657 |
| h=.00625 matched defect | 1.124 | 3.603 | 9.059 | 18.682 |
| h=.003125 matched defect | 1.124 | 3.607 | 9.146 | 20.069 |
| h=.003125 continuum defect | 1.124 | 3.616 | 9.454 | 32.590 |
| h=.003125 matched vacuum | .254 | .574 | 1.210 | 2.332 |

The coarse continuum 64 entry is outside the original h*Omega<=.5 bound. At finest matched lattice, M0=.021105 by 64 and M2=.105293, vs .099933 at 32. M4 more than doubles from 32 to 64; P1 also grows (11.725 to 25.658). Massless P2 matched defect M4=10.556 to 22.999. Matched vacuum grows more weakly. At Omega=64 there is still meaningful grid dependence even in the matched diagnostic (P2 18.682 to 20.069); do not assert an exact continuum UV divergence or asymptotic power. Nonetheless the tail is not solely the coarse numerical artifact, and saturation of unweighted response does not guarantee saturation of frequency-weighted response.

## Interpretation / next step

Finite spatial size gives a natural *elastic* resolution scale omega~c/R for a massless incident field, but rapidly oscillating incident waves can access rapidly oscillating internal modes; small coherent form factor does not bound the full excitation spectrum. A Gaussian also rolls off, hence this is not yet topology-specific. A physically specified field and coupling, full spectral response (angular channels, absorption/scattering and spectral shifts), and its actual vacuum frequency measure are needed before asserting a finite ZPF dressing. Neither radial overlaps nor finite particle size cut off vacuum-only loops responsible for bulk induced gravity. Keep that question separate from the GC30 common-metric problem. No physical frequency in Hz without calibrated units and propagation speed.

## Reproducibility / budget

Run `python src/exgc29a.py` then `python src/exgc29a_dispersion.py`. Raw compact results massless.json, massive.json, dispersion.json. Main production 13.11 CPU s; first attempt completed calculations but failed JSON serialization of numpy bool, corrected conversion only and reran (~13 s initial); diagnostic 25.98 CPU s; total approximately 52 CPU s (<.02 core-hour), well inside 2 core-hour cap. No changed scientific criteria.
