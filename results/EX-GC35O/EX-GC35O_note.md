# EX-GC35O — reciprocal local clocks with wave and carrier sources

2026-11-24 | EXPLORATORY | ChatGPT | [Plan](../../spec/EX-GC35O_plan.md). Started after GC35L and GC35M reported. PI closure pending.

**Two-way feedback worked in the reduced event host. The competing token histories remained inequivalent.** This is a scalar static-response reduction of Y1 on a periodic line. GC35M did not qualify a full local tensor evolution, so the toy imports Y1's static lapse symbol and the isotropic spatial relation, with a separate assumed retarded lapse equation.

## Local rule and controls

Evolve a packet q and its momentum p with q̇=(N/A)p and ṗ=div((N/A)grad q) in one spatial dimension. The retarded perturbation n obeys n̈=Δₕn−s(ρ−ρ̄initial) with s=.18 and unit propagation speed. Its nonzero-mode static symbol is n̂=−sρ̂/|K|², with the same spatial symbol as static Y1. Set N=1+n and A²=1−2n in the weak regime. The packet's local energy density sources n and the clock rate at a site sets its next coordinate-time interval Δτ/N. The zero-mode background subtraction is fixed by the initial data. The retarded evolution and the algebraic A are imposed in this toy and do not reproduce the dynamical Y1 tensor constraints. The static sinusoid control had a relative Poisson residual of 2.55e−14.

Use 48 sites, spacing .25 and a coordinate interval through t=6. Each event reads prior versions and linearly predicts a neighboring field at its current timestamp. A neighbor's published next deadline certifies that all versions at or before that timestamp are available. The discrete integrator is a first-order explicit update in the dynamical variables. No globally synchronous step or spatial projection is applied. The reference executes events in chronological order and randomized ready schedules may advance some neighborhoods while others lag.

| Proper step | Field events | Per-site counts | End lapse range | Matching randomized schedules |
|---:|---:|---:|---:|---:|
| .04 | 7,190 | 144–155 | .9291–1.0895 | 8/8 |
| .02 | 14,409 | 289–311 | .9326–1.0902 | 8/8 |

Every match is byte-for-byte across fields, timestamps and version dependencies. No lapse clamps were reached. At the finer timestep, disabling the source removed the lapse perturbation, and disabling feedback while still evolving the sourced lapse changed the final packet by 0.4467 in the same unscaled spatial L² norm. Freezing the lapse at the initial static solution changed the packet by 0.2586. Halving the dynamic timestep changed the packet by 0.0620 and lapse by 0.0177. These are finite-window sensitivity checks, with no precision convergence gate.

## Contested token fixture

Embed the GC35L paired redistribution on three marked sites and two labelled sinks. The local token counts add to the lapse source, and the current lapse sets the next proper-time firing deadline for an eligible site. For the hybrid fixture the executor uses chronological order, deterministic site order for equal timestamps, and two distinct legal policies for selecting a pair of tokens. In particular the center can use its original tokens or consume an arrival from the left. The paired rule conserves seven total tokens including sinks. It ends at (2,1,1,1,2) for both policies, with identical fields and hop timestamps but inequivalent labelled carrier DAGs. The last hop occurred at t=1.36793 with feedback and t=1.36000 with the clocks held at unity. Adding carriers to the packet source changed the final lapse by 0.2960 in the same spatial L² norm.

These runs couple the clock and the sources in both directions. Schedule equality of the versioned field implementation and invariance of graphs under contested token selection are separate observables. The latter failed even though the matter ledger and the metric agreed for this pair of branches. Full tensor Y1 reciprocity, a joint exact momentum law and broad competing-event invariance remain open. [Source](../../src/exgc35o.py) · [Raw results](results.json)
