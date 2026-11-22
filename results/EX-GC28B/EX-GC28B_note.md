# EX-GC28B — connected vacuum interaction (EXPLORATORY)

**Summary:** Two disjoint positive compact scatterers in a free real scalar vacuum attract through connected quantum fluctuation energy, without a linear classical source. In 3D the potential decays as d^-3, hence force d^-4, and its amplitude is controlled by scattering lengths rather than rest energies. In 1D positive barriers approach the familiar -pi/(24d) energy irrespective of barrier strength. A small mediator mass exponentially suppresses both. This is a successful Casimir control, negative for Newtonian gravitation from this particular quadratic vacuum coupling. It does not rule out induced metric dynamics, a collective gapless pole, topology-dependent couplings, or other quantum vacuum mechanisms.

Design, implementation, run: ChatGPT. PI authorized; independent review and PI closure pending. Pre-run plan `spec/EX-GC28B_plan.md`, infrared correction `spec/EX-GC28B_quadrature_check.md`, review disposition/Born derivatives `review/EX-GC27AB_28A_response.md`. 135 initial rows + 55 infrared rows; 36.77 CPU-s (<.011 core-h), 1 thread, within 2 core-h. Reproduce with `python src/exgc28b.py` then `python src/exgc28b.py infrared`. The compact records contain every computed energy and parameter.

## Method and controls

H=-Laplacian+m²+V1+V2, V positive and supported in unit ball, compact truncated Gaussian sigma=.35 or .5 and integral Q. No linear source or direct potential between defects. At disjoint support, subtract each isolated defect and add vacuum, so connected energy is UV finite; no subtraction of interaction terms. Use free-space Euclidean Green function and Cholesky whitening of each defect self-block, integrate log det(I-BB^T)/(2pi) over imaginary frequency. Equivalent to the scalar TGTG approach (Kenneth & Klich, arXiv:0707.4017; Emig et al PRD 77 025005). 3D diagonal Green value is equal-volume sphere cell average; finite spacing is an approximation, controlled by h refinement. No artificial enclosing box. Positive V ensures stability. The direct finite lattice eigenfrequency sum independently matches its own determinant integral to 6.44e-14; point delta exact test to 2.2e-19; zero scatterer exactly zero. These checks validate the connected subtraction, not the physical model.

## 3D massless results

| d | Q=.05 | Q=1 | Q=10 |
|---:|---:|---:|---:|
| 4 | -2.02831e-8 | -6.53713e-6 | -1.79121e-4 |
| 8 | -2.45436e-9 | -7.82828e-7 | -2.02404e-5 |
| 16 | -3.04371e-10 | -9.66329e-8 | -2.44199e-6 |
| 32 | -3.79693e-11 | -1.20283e-8 | -3.00922e-7 |
| 64 | -4.74362e-12 | -1.50114e-9 | -3.73790e-8 |

h=.25, sigma=.35. Effective potential exponents over 32→64: 3.00077, 3.00231, 3.00909 respectively. Scattering lengths .00395259, .0702765, .349898. The energy divided by flight-distance asymptotic -a1*a2/(4pi*d³) at d64 is 1.00022, 1.00127, 1.00576. At h=.2 changes vs h=.25 at d64 are -.00616%, -.1067%, -.4734%; quadrature 48→96 and cutoff 32→48 changes selected 3D values by <2e-12 relative. Mixed Q1=1,Q2=10 gives -7.49069e-9 at d64, close to geometric mean of equal-species values (7.49057e-9). Changing sigma of one Q=1 scatterer .35→.5 changes d64 energy from -1.50114e-9 to -1.52984e-9 (~1.9%) at fixed Q. Thus factorization into scattering lengths is supported, not proportionality to inertial energy (not measured). At Q=10 the amplitude is not proportional to Q²: compared with Q=1 an increase of 10x both strengths raises attraction by ~24.9x, not 100x. The weak Born coefficient is recovered only at weak coupling. The potential d^-3 gives force d^-4; not Newtonian. When m=.15, Q=1 at d64 is -3.97682e-17, only 2.65e-8 of the massless energy; consistent with exp(-2md) times power prefactor.

## 1D infrared check

Initial massless quadrature did not converge near the logarithmic infrared endpoint. Preserved initial records; prospectively registered t=u² substitution and 96 vs 192 nodes. For Q=.05 at d4 the original value -0.00382329 changes to -0.00416492; 192-node -0.00416490. All selected refined relative differences <5.1e-6. Grid h=.1 vs .05 changes chosen cases <.2%. Refined massless exponents and amplitude:

| Q | p at 32→64 | p at 512→1024 | E(d1024)/[-pi/(24d1024)] |
|---:|---:|---:|---:|
| .05 | .605992 | .949113 | .963233 |
| 1 | .967181 | .997824 | .998489 |
| 10 | 1.01269 | 1.00078 | 1.00054 |

Weak 1D potentials are not uniformly weak in the IR; any positive barrier ultimately reflects long wavelengths, giving -pi/24d irrespective of Q. The apparent inverse-square force in 1D does not produce Newtonian gravity in 3D and its asymptotic coefficient does not track particle mass. Massive controls decay exponentially. No comparison with real inertial mass; this is a static scattering model with inserted potentials rather than dynamical topological defects.

## Programme implication and proposed next direction

Separate a **direct** two-fluctuation exchange from an **induced** propagating collective field. For quadratic local scattering the 3D connected loop contains two massless propagators and a frequency integral, producing d^-3 rather than single massless static Green's function d^-1. Sakharov induced gravity is not a claim that ordinary Casimir forces equal gravity: it proposes vacuum contributions to a metric action. A useful next screen is whether a dynamical vacuum-order parameter has a stable gapless mode with a static 1/k² response and whether defect energy sources it through a common action. Derive the dispersion, sign and stiffness, source coefficient including quantum/counterterm contributions; compare with 26C/F which assumed a massless Nordström scalar. Do not introduce a massless mediator and then claim its existence emerged. If an odd/linear source is inserted, treat its 1/r interaction as an explicit control. Quantum entanglement can be investigated as a correlate of the connected vacuum state, not assumed to supply an additional force. A genuine 3D 1/r interaction and universal coupling to total energy remain open.
