# EX-GC32C — Vacuum-selected relational geometry

2026-11-22 | ChatGPT | EXPLORATORY | complete; independent review and PI closure pending. **PI CLOSED 2026-11-22** (Richard: "I sign off on GC31A–GC33A."). Claude review: review/EX-GC31A-33A_claude_review.md.
[Registered plan](../../spec/EX-GC32C_plan.md), [Python source](../../src/exgc32c_relational.py), [exact/reference/perturb data](geometry_reference.json), [search data 1](geometry_search_1.json), [search data 2](geometry_search_2.json), [excitation data](excitation_response.json). Reproduce with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python src/exgc32c_relational.py. 15.03 CPU seconds and 15.04 wall seconds; single thread, well inside 2 core-hour cap. All registered runs completed, no parent changes and no further production. The plan and LOG were committed before production.

## 1. Disposition

**NOT SUPPORTED:** ordinary positive harmonic zero-point energy as the selector of extended homogeneous low-dimensional relational geometry. Both an exact monotonicity obstruction and a controlled fixed-budget clustering preference are found. Introducing conserved oscillator excitation energy modifies the link-gradient in a state-dependent manner but does not create a universal total-energy source. This is a restricted model result, not a no-go for emergent spacetime or Sakharov induced stiffness. A geometry-selecting energy is not the same as gravity on a specified metric.

## 2. Model and analytical result

All unordered abstract label pairs may become links. No coordinates, lattice distances or target geometry enter the energy or swap proposal. For graph G:
    K(G)=mu² I + L(G)
    H_G=(p^T p + q^T K q)/2
    U0(G)=tr sqrt(K)/2
mu>0 is a supplied gap/regulator, link stiffness unit. Link variables are classical slow parameters here (Born–Oppenheimer potential), not quantum dynamical degrees of freedom. The search algorithm is an optimization/control, not physical time or an energy-conserving joint link/oscillator dynamics. N, clock, oscillator normalization, and edge constraints remain inputs. A bare term proportional to number of edges cannot distinguish equal-edge controls.

A link between i,j contributes w vv^T to K, v=e_i-e_j. For positive K, differentiating under the trace gives
    dU0/dw = (1/4) v^T K^(-1/2) v > 0.
Thus every nonzero link strictly raises U0 relative to deleting it. The empty graph is the unique unconstrained minimizer (at fixed N, mu>0). A positive sum over independent oscillator sectors of this form retains the obstruction. A regulator modification preserving this monotonic structure does not solve it. Trace sqrt is concave on positive matrices; at fixed total Laplacian trace the concavity tends to favor greater spectral spread, not homogeneous eigenvalues. This doesn't order every pair of graphs (would require majorization, and degree/triangle effects also matter).

For large mu relative to maximum sqrt(lambda_L), expansion:
    U0 = N mu/2 + tr L/(4mu) - tr L²/(16mu³) + tr L³/(32mu⁵) + ...
At fixed edge count tr L=2m and tr L²=sum d_i²+2m. The second correction favors unequal degrees. At fixed regular degree d, tr L² fixed; tr L³=Nd³+3Nd²-6T (T triangles), so the first distinguishing term favors triangles; higher moments determine remaining differences. This explains why simple regularity is not enough. The tested small-gap behavior is calculated directly, not inferred from this asymptotic expansion.

## 3. Exact graph population
Every bitmask among 2^15=32,768 labeled simple graphs at N=6 included. 93 class summaries (31 per regulator). Empty unique global minimum verified at mu=.1,1,5. At fixed edge count, 11/16,11/16,10/16 unconstrained class minima disconnected respectively (includes empty graph; dense edge counts eventually force connectivity). Some graph minima change with mu. Connectedness required by hand does not guarantee evenly distributed geometry.

Selected exact mu=.1 minima (energies for whole graph):
| Class | Unrestricted minimum geometry | U0 | Connected minimum geometry | U0 |
|---|---|---:|---|---:|
| 5 edges | K4 minus edge, 2 isolated | 2.861371 | star K1,5 | 3.285740 |
| 6 edges | K4, 2 isolated | 3.153748 | degree pattern 1,1,1,2,2,5 | 3.650714 |
| 9 edges | K5 minus edge, 1 isolated | 4.324922 | degree pattern 1,2,3,3,4,5 | 4.607532 |
| 12 edges | degree pattern 2,4,4,4,5,5 | 5.448705 | same | 5.448705 |
| 2-regular | 2 disconnected triangles | 3.569870 | C6 | 3.791172 |
| 3-regular | triangular prism | 4.732110 | triangular prism | 4.732110 |

Counts of tied labeled graph minima and representative invariant/spectral types retained in archive (tolerance 1e-10); e.g. reg2 unrestricted 10 labeled partitions, connected 60 cycles, reg3 60 prisms. Tied-type examples are grouped by spectrum/degree/triangle/component signature, not a general isomorphism classifier. At mu=5, 5-edge minimum changes to star and connected 9-edge minimum to degree pattern 2,2,2,2,5,5; disconnected triangles beat cycle in all three regulators. These are exhaustive small-N minima only; not evidence of an infinite-volume phase.

## 4. Fixed-degree large graph controls and search
Four size families, degree four (N=20,40,60) compared with periodic 4x5,5x8,6x10 reference, degree six N=64 compared with periodic cube 4x4x4. Disconnected dense blocks K5 or K8 minus a perfect matching; connected necklace (one edge per block removed and endpoints joined cyclically), two randomized regular connected controls independent seeds. All degrees verified, all spectra invariant under an independent relabeling to 1e-10. These controls have the same edge budget and degree pattern; no energy term distinguishes them except the computed spectrum.

mu=.1, U0/N:
| N, degree | Extended lattice | Connected necklace | Disconnected blocks | Random17 | Random61 |
|---|---:|---:|---:|---:|---:|
| 20, 4 | .953283 | .934230 | .905321 | .952130 | .952175 |
| 40, 4 | .957455 | .934996 | .905321 | .958072 | .958334 |
| 60, 4 | .958547 | .935061 | .905321 | .959925 | .960347 |
| 64, 6 | 1.192001 | 1.165814 | 1.149794 | 1.192048 | 1.191788 |

All 12 family/regulator combinations satisfy U_blocks < U_necklace < U_lattice and both random controls above U_necklace; even mu=5 where differences small but resolved. This directly fails the registered positive preference gate, independent of optimization quality. The torus can beat random controls sometimes, but that is inadequate. Necklace global clustering .7, lattice 0, blocks 1 (.8 in degree6). Degree4 necklace diameter 6,12,18 with gap .26795,.07237,.03258, versus lattice diameters 4,6,8 with gap 1.38197,.58579,.38197. A large graph diameter or apparent one-dimensional slow mode alone would mistake a chain of cliques for homogeneous geometry.

48 bounded search runs, 800 valid swaps each (38,400), from lattice, random17 and necklace, connectivity enforced or allowed to fail, mu=.1/1. All 32 lattice/random starts improved their best energy; all 16 necklace starts retained their initial best (not a proof of global or local optimality). Best clustering >.1 in all 48 runs; mean .497 and .507 connected-only at mu=.1/1, .521 and .501 unrestrained. Only 1/24 best graphs in unrestrained search fragmented (mu=.1, from random); most remain connected, so do NOT interpret finite search as proving dynamical disconnection. Disconnected block comparator is lower independent of search; barriers and limited swaps matter. 8 independent lattice perturbation assays (1024 valid connected proposals) yielded 26 to 94 descending moves out of 128 (20.3%–73.4%), each exceeding registered <10% threshold. Neither selected homogeneous geometry nor metastability established for the toroidal references in this move set.

## 5. Excitation forces from same energy functional
For a fixed state of one normal-mode quantum, independently track source excess energy omega_j=sqrt(lambda_j(K)). In the unperturbed eigenstate, fixed-state link variation contributes
    dE_src/dw = (u_j(i)-u_j(j))²/(2 omega_j).
This is the expectation variation in the chosen state, including within a degenerate eigenspace; not an asserted uniquely tracked eigenbranch at a crossing. Source force on w is **minus** this derivative, and reduces the preference for stiffness along the excited differences; no derived force on a particle, no gravity or inverse-square inference. Every possible pair including absent links is eligible for this variation (centered finite differences only on present positive links).

36 graph/regulator cases, 108 source mode cases, 72 coherent controls. FD vacuum derivative max relative error 3.63e-8; fixed-state source expectation derivative max abs error 5.56e-11 (relative error for zero derivative undefined); uniform source gradient max 1.06e-27 (active and absent); coherent ledger max abs error 1.56e-15. All analytical-control criteria pass. At each regulator, uniform mode carries positive energy mu, yet has exactly zero extra link force because all coordinates oscillate together. This provides a structural witness against **all energy sourcing link deformation** in this fixed-gap scalar model. It is a delocalized excitation, not a massive particle; for mu→0 its single quantum energy vanishes, so do not make this a statement about a massless regulator-free theory.

Nonuniform modes generate gradients. For the N40d4 lattice at mu=.1:
| Mode | E_src | Mean active-edge derivative increment | Mean absent-pair derivative increment |
|---|---:|---:|---:|
| uniform | .1 | ~0 | ~0 |
| first nonzero | .77187 | .004743 | .036474 |
| highest | 2.76189 | .017239 | .008375 |

These depend on the mode shape and normalization, not merely source energy. Discrete rank comparison uses surfaces U0+omega_j with fixed index per graph (not source adiabatic continuation). Necklace remains lowest among the three connected references for all 12 families/regulators and all three occupation choices. Low mode reorders lattice/random in one group, high mode in two; does not rescue geometry. The independently computed source energy is included consistently, not inferred by fitting link forces.

Matched-localization coherent states: input q displacement at vertex zero vs input p impulse there, each E_src=1 in the same H. Excess link derivative in the coherent state is (q_i(t)-q_j(t))²/2. For N40d4 lattice mu=.1 the initial maximum active-edge increment .249377 for q versus zero for p; over fixed t=0..40, time-average mean active increment .006246 vs .006064 and absent increment .006242 vs .009371. Both conserve energy (all cases). Phase/preparation matters; time averages need not be equal because local q vs p prepare different mode-energy distributions. This is not a test of relativistic equivalence or proof that a stress-based source cannot exist. It shows a scalar density/strain link force is not total-energy universality. The nonlocal absent-pair forces are additionally problematic for physical locality unless link dynamics supplies it. No link evolution or source localization survival shown.

## 6. Big picture and recommended disposition

This result connects to GC32A/B: response, ground-state data, and total energy are distinct information. A positive scalar oscillator vacuum energy is nonzero but not sufficient to select physical spatial geometry; imposing a fixed valence and connectivity does not solve the problem. GC28D's Sakharov induced curvature stiffness on supplied geometry remains logically independent. To get a gravity mechanism we need a stable relational vacuum with emergent locality, robust common response and energy-sensitive geometry sourced by the complete conserved stress-energy, not just link-strain expectation. A further search that merely adjusts the scalar spectral energy or hand-selects short-cycle targets would risk repeating GC31.

Prior art is highly relevant: Konopka, Markopoulou and Severini, [Quantum Graphity: a model of emergent locality](https://arxiv.org/abs/0801.0861), Phys Rev D 77,104029 (2008), use quantum graph variables and competing valence/loop terms to investigate low-energy ordered geometry. Hamma et al., [A quantum Bose-Hubbard model with evolving graph as toy model for emergent spacetime](https://arxiv.org/abs/0911.5075), explicitly make interaction structure a quantum dynamical variable. GC32C does not reproduce these models or demonstrate that they fail: the tested vacuum functional is much narrower.

Recommendation: stop pure-positive harmonic ground-energy minimization rather than increase N or tune an assumed lattice. Next consider **GC33A (paper screen only, proposed)**: can an independently motivated vacuum *dynamical phase* select spatial locality and support a common mode response? Inspect quantum graphity / quantum matter-link exchange and genuinely quantum collective constraints/entanglement, including competing vacuum spectral contributions with independently constrained signs, rather than fitting valence/cycle coefficients to geometry. Ask what is conserved in a joint link/matter Hamiltonian, whether locality/dimension and common cone survive perturbations, and what energy/stress drives changes. Require an analytical selection mechanism or discriminating negative prediction before simulation. This is a change of mechanism, not more precision on GC32C. Avoid claiming a new theory from the graphity prior art.

**State:** all registered numerical assays done, selection hypothesis NOT SUPPORTED; analytical controls passed. Pending independent review/PI closure. No GC33A or subsequent production authorized.
