# EX-GC27A/B — static vacuum-scale sensitivity

ChatGPT, 2026-11-22. PI approved GC27. Reviewed (review/EX-GC27AB_28A_claude_review.md; response review/EX-GC27AB_28A_response.md). **PI CLOSED 2026-11-22 (Richard: "go ahead and close 27A/B")**. Registered plan `spec/EX-GC27AB_plan.md` and LOG posted before production. All registered static cases completed in 1444.13 CPU-s = 0.4012 core-h (<2). No new gravitational field, force, or quantum fluctuation dynamics was run. Skyrme reference model, not Model Gc.

## Main result

The scale identity explains universal scalar sensitivity in the massless model; a second independent scale removes that guarantee. The specific single and bound-pair states tested have very similar energy partitions, so the resulting species spread is small and diminishes over the two grids. This does not establish a new universal principle, nor a continuum violation. It provides an explicit way to diagnose what must be explained by a microscopic vacuum theory.

Let E=A E2+B E4+C V0, eta=ln(chi/chi0), exponents (a,b,c), normalized derivative k=(a E2+b E4+c E_V)/E at eta=0. Here E2,E4,E_V include their baseline coefficients. Write v=(E2-E4+3 E_V)/E and p=E_V/E. An exact algebraic identity for the recorded values is

`k = (a+b)/2 + (a-b)*v/2 + (c-2*a+b)*p`.

The middle term is the measured virial residual; removing it is an algebraic diagnostic, not a continuum extrapolation. For a stationary continuum state v=0. The remaining sensitivity is universally fixed only if the scale weights are compatible, or if relevant energy fractions happen to agree.

## A: massless control

| h | exponents (a,b) | k B1 | k B2 | continuum |
|---|---|---:|---:|---:|
| .2 | 2,0 | .935938 | .948418 | 1 |
| .2 | 0,2 | 1.064062 | 1.051582 | 1 |
| .2 | 2,-2 | -.128125 | -.103163 | 0 |
| .2 | 2,2 | 2 | 2 | 2 |
| .15 | 2,0 | .963826 | .970375 | 1 |
| .15 | 0,2 | 1.036174 | 1.029625 | 1 |
| .15 | 2,-2 | -.072347 | -.059249 | 0 |
| .15 | 2,2 | 2 | 2 | 2 |

The full species difference is the measured virial term by identity. The size-only control changes particle scale without sourcing in the continuum. The raw binding sensitivity for quadratic scaling is .630/.812, compared with 1 in the continuum; subtraction amplifies the finite-lattice errors. Common scaling gives binding sensitivity exactly 2. Parity reverses degree and changes energy by <1.5e-13. Vacuum is exactly zero energy and projected gradient.

## B: independent positive potential scale

The coefficient was fixed once at C0=.17839644285107545 (m=.29866071289, from registered mR=.3 using coarse B1 radius 1.00448431), before B baselines. No species-specific tuning. Four relaxed states retained degrees -1/-2. The potential shrinks B1's coarse radius from 1.00448 to .92099 (8.31%); binding fraction remains about 4%.

| h | case | k B1 | k B2 | B2−B1 raw | B2−B1 excluding virial term |
|---|---|---:|---:|---:|---:|
| .2 | scale-consistent (2,0,4) | .931703 | .944834 | .013131 | 0 |
| .2 | independent (2,0,0) | .863539 | .874036 | .010496 | -.002635 |
| .2 | potential-only (0,0,1) | .017041 | .017700 | .000659 | .000659 |
| .15 | scale-consistent (2,0,4) | .964544 | .970477 | .005933 | 0 |
| .15 | independent (2,0,0) | .890613 | .894771 | .004158 | -.001775 |
| .15 | potential-only (0,0,1) | .018483 | .018926 | .000444 | .000444 |

Independent scale, after virial separation: B1/B2 .931836/.929201 at h=.2 and .926069/.924294 at h=.15. The common ~7% source suppression is much larger than the between-species difference. A common suppression alone could be absorbed into an overall coupling normalization, so it is not a violation of universality. The potential fraction differs by only .000659/.000444 (absolute), and that difference shrinks with refinement. Accordingly the continuum composition dependence of these two configurations is unresolved. The potential-only source has near-zero binding response despite nonzero binding energy, a useful contrasting control. The raw independent binding sensitivity is .603/.793; not a continuum result. Full data in `summary.json`.

## Checks and provenance

- A inputs: B1 arrays reused from hash-verified B2 local states (not the exact E1/E2 arrays); coarse B2 regenerated bit-identical to E1 hash 04ffce6a..., fine B2 regenerated on Linux with energy difference -1.76e-10. All four total energies agree with E1/E2 well within .1%; fine B1 E2 component differs ~.100% though total difference ~1.36e-6 relative. Prior and actual state identity distinguished.
- Weighted gradient validation relative error 1.77e-9; vacuum exact. Maximum fixed-shape central derivative minus analytic k .000534 (expected finite perturbation error). Maximum relaxed central derivative difference .000197, across twelve coarse continuations. Every relaxed endpoint degree retained; largest projected gradient 1.07e-5. Baselines B max gradient 1.20e-6. Optimizer iteration-limit exits recorded, not hidden. Outer 1.5-unit energy fractions <.15% (A), <.054% (B); no box convergence claim.
- Fixed and rescaled shape sweeps at eta ±.02 and ±.05, A; B ±.02. No re-tuning or failed topological states discarded. All baseline hashes and raw numerical records included; arrays local and reproducible from code. `python src/exgc27_summary.py` independently recombines recorded source/virial/binding data.

## Recommendation

Do not spend the next phase chasing a 0.2% composition difference in these near-similar partitions. The important result is the exact scale condition. A truly universal source requires a common scaling principle across all independent contributions; topology alone gives no such guarantee. A larger topologically distinct cohort may eventually test this, but first pursue Candidate 4: whether vacuum fluctuations supply or constrain such a scaling principle. Start with a small quantized harmonic normal-mode calculation with a known reference and controllable cutoff, then the knot Hessian/kinetic metric. Measure vacuum-subtracted energy sensitivity and connected pair interaction; no assumption that either is gravitational. GC27C (dynamical vacuum response) can follow a motivated vacuum action rather than merely reintroducing Nordström. PI decides.
