# EX-GC30B — Common effective geometry robustness (EXPLORATORY)

**PI CLOSED 2026-11-22 (Richard: "I think those runs can be closed"; Claude review review/EX-GC30A-B3_claude_review.md).**

Two coupled condensates support one gapped and one gapless excitation sharing the same effective low-energy dispersion geometry at a tuned reference vacuum. Both respond identically to shared anisotropy and flow at that density. But a common density change alone splits their spatial metric coefficients, without breaking exchange symmetry; only 4 of 15 registered symmetry-preserving test conditions pass the 1% gate. This is a negative universality result for this analogue benchmark, not a no-go for all emergent metrics. Recommend moving to collective tetrad models rather than GC30C energy sourcing on this failed candidate.

Design/run: ChatGPT, PI authorized 2026-11-22; exploratory, preregistered spec/EX-GC30B_plan.md (e13555f4), LOG before production (4838809a). No closure sign-off. Reference BEC, not Model Gc. Two core-hour cap, one thread; ~2.33 CPU s incl. invalid first attempt. No subsequent stage run.

## Model and validation

Hamiltonian density sum_a [(grad psi_a)* M_a (grad psi_a)/2+V_a n_a+g_a n_a^2/2]+g12 n1 n2−lambda(psi1*psi2+c.c.); hbar=1. Background sqrt(n_a)exp(i p.x−i mu t), mu=0, V_a chosen to satisfy stationary GP equation (not fitting propagation). Positive densities/inverse masses, stationary residual <1e-12. Linearization: eps_a=k^T M_a k/2, s_a=p^T M_a k, A_aa=eps_a+g_a n_a+lambda sqrt(n_b/n_a), A_ab=g12 sqrt(n_a n_b)−lambda, B_aa=g_a n_a, B_ab=g12 sqrt(n_a n_b). BdG matrix [[A+S,B],[-B,-A+S]], S=diag(s). Positive-frequency roots must have positive Bogoliubov norms; all sampled 23 cases pass, not a proof of global stability. No curved metric supplied. For equal components the eigenvectors separate into in-phase and out-of-phase sectors; at lambda>0 the gapless and gapped branches are distinct. Sorted roots sufficient in low-k window, no species identification in degenerate decoupled control. Relative-case fits describe sorted branches, not independent microscopic species labels.

At common density n, g1=g2=g, M1=M2=M, exact squared comoving frequencies are eps(eps+2n(g+g12)) and (eps+2lambda)(eps+2lambda+2n(g−g12)), eps=k^T M k/2. This follows by diagonalizing in symmetric/antisymmetric field basis. Therefore the quadratic matrices are B_+=n(g+g12)M and B_-=[n(g−g12)+2lambda]M; gap²_-=4lambda[lambda+n(g−g12)]. Common geometry at quadratic order iff lambda=n*g12. Neither exchange symmetry nor positivity enforces that equality. See Liberati, Visser, Weinfurtner (2006), https://arxiv.org/abs/gr-qc/0510125 and Sindoni (2010), https://arxiv.org/abs/1011.4411. These demonstrate tuned/symmetry-limited analogue geometry; our numerical test checks robustness with specified controls, not new physics beyond those works.

Fit even frequency squared to gap² + k² directional coefficient + k⁴ term; odd frequency to linear and cubic terms, signed k in 10 directions (9 fit, 1 holdout), 8 samples per sign per direction, two nested windows .01–.08 and .005–.04. Reconstruct full symmetric spatial matrices and drift vectors independently from spectral roots. Mismatch D_B=||B_+−B_-||_F/mean(||B_a||_F); drift mismatch normalized to mean characteristic speed, gate 1% both. Max window uncertainty 1.08e-7 (B), 1.25e-8 (drift); maximum holdout errors 4.64e-9 and 8.82e-10 respectively; all below .002. Symmetric controls match exact omega² within 4e-15 and B within 1e-10 relative. All corrected spectra real within tolerance, two positive roots and positive Bogoliubov norms. This is homogeneous/local dispersion analysis, not a simulation of spatial curvature or moving particles. A gapped branch's low-k curvature is not its exact high-frequency characteristic cone.

## Results

| Condition | D_B (%) | Notes |
|---|---:|---|
| Decoupled equal | 0 | two gapless positive control |
| Exchange interaction .1, lambda=0 | 20 | exact negative control |
| Baseline lambda=g12=.2 | ~0 | gap²=.8; B=1.2 I |
| n=.8 | 8.000 | fixed lambda |
| n=.9 | 3.636 | fixed lambda |
| n=.95 | 1.739 | fixed lambda |
| n=1.05 | 1.600 | fixed lambda |
| n=1.1 | 3.077 | fixed lambda |
| n=1.2 | 5.714 | fixed lambda |
| g12=.18 / .22 | 3.333 / 3.333 | fixed lambda, n=1 |
| lambda=.18 / .22 | 3.390 / 3.279 | fixed g12, n=1 |
| Common flow | ~0 | fixed n=1 |
| Common anisotropy | ~0 | fixed n=1 |
| Common flow+anisotropy | ~0 | fixed n=1 |
| Common flow+anisotropy, n=1.1 | 3.077 | fixed lambda |
| Relative density | .102 | deliberately broken exchange |
| Relative g | 1.667 | deliberately broken exchange |
| Relative tensor | .00000138 | deliberately broken exchange |
| Relative tensor+flow | .119 | deliberately broken exchange |
| Retuned n=.8 / 1.2 | ~0 / ~0 | lambda=.2n; excluded from gate |

All drift mismatches <3e-9 in fine fits. 4/15 symmetry-preserving fixed/parameter test cases pass, 11 fail (6 density, 4 parameter, 1 mixed density/flow/anisotropy); controls, relative cases and retuned cases excluded from gate. Common anisotropy gives same B entries: diagonal (1.38,1.02,1.32), xy=.144, yz=.06; xz≈0. Shared anisotropy and flow at tuned density are genuine conditional kinematic successes. Relative tensor leading equality is not evidence of symmetry protection: at tuned point the first k² perturbation samples the average kinetic tensors in the symmetric/antisymmetric eigenstates, with mixing at higher order. Other relative controls are not necessarily large; relative density gives small mismatch; do not overstate negative results.

At n=1 fixed couplings, gapless c²=1.2, gapped c²=1.2. At n=1.1, c²=1.32 and 1.28; mismatch 3.077%. Both experience identical physical density perturbation yet differ. Exact formula D_B=2|n*g12−lambda|/(ng+lambda) for shared M, since the common tensor norm cancels. At n=.8 and 1.2 retuning lambda=n*g12 restores equality exactly but changes microscopic Hamiltonian; not a universal response. A hypothetical dynamical law enforcing lambda=n*g12 is a new hypothesis, not a finding here; must work for arbitrary background gradients and independent probe sectors, with symmetry/radiative stability demonstrated.

High k quantum pressure: baseline along x at k=2 full frequencies 2.96648 and 3.09839 versus quadratic approximation 2.19089 and 2.36643. A common quadratic metric is approximate even at the calibrated point. Nonrelativistic microscopic rest frame remains, no tensor gravitational pole, no universal nonlinear stress sourcing derived. C4 full metric not assessed; gaps vary with n, so common propagation speed does not by itself establish clock normalization. No soliton dynamics or gravitational sourcing measured.

## Decision

Preregistered universality gate fails. Do not proceed to GC30C on this reference. Recommend a separate GC30B2 or GC30T analytical/numerical test of emergent tetrads around distinct Weyl nodes or inequivalent excitations in a specified common microscopic Hamiltonian. Require shared effective inverse metric under deformation without assigning a common metric, and distinguish symmetry-related nodes from genuinely different modes. A microscopic common hopping rule alone does not demonstrate universality between bosonic and fermionic sectors. Only after a robust kinematic candidate passes, test energy-induced backreaction and massless tensor dynamics; preserve Sakharov curvature stiffness as downstream mechanism. No additional numerical precision needed on this exact counterexample.

## Validation deviation and reproducibility

First attempt 1.11 CPU s invalid: integer B array discarded off-diagonal interaction; exchange_split complex spectrum; baseline exact omega² error .44. Logged prospectively before rerun in spec/EX-GC30B_validation_correction.md (535eedfa). Explicit float n,g,M,p and enforcing analytic control assertion are the only changes. First-run records retained in local checkout results/EX-GC30B_failed_dtype_attempt, diagnostics in correction note; they are not scientific evidence. Corrected run 1.21 CPU s. Source src/exgc30b.py, run `python src/exgc30b.py`; compact case records cases_0.json through cases_4.json contain raw directional spectra, fitted matrices, gaps, quartic terms, holdout, stability and higher-k samples. Same-author analytic recomputation is not separate-author replication. Review and PI closure pending.
