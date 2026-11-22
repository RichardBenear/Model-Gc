# EX-GC28F — shielding and anti-shielding (EXPLORATORY)

**Summary:** A shell that changes only the secondary-field stiffness changes the local field but cannot shield a fixed integrated monopole: Gauss flux is conserved. A polarizable shell that adds an explicit local phi² term strongly screens (or enhances for mild negative response) the far field, even while the exterior remains massless and Coulombic. With self-consistent core feedback a stiffness shell can also indirectly change the monopole by modifying the core response. This is not GR gravitational shielding; the modeled material response explicitly couples to a particular secondary scalar, not total energy.

Design/run ChatGPT, review and PI closure pending. Plan `spec/EX-GC28DEF_plan.md`; source `python src/exgc28def.py F`. 44 rows, all converged, single thread; final run ~6.42 CPU-s plus initial similar run before fixed-profile radial eigen diagnostic was included (total ~13 CPU-s, well below 2 core-h). Baseline C=.178396, q=.25, alpha=1, Q=.0108263, phi(0)=.140662. Smooth shell profiles sin², radii 2–4 and 4–8, outer R24/48. Fixed source is the *net baseline* q rho-alpha rho phi, not just bare q rho, keeping it unchanged while changing shell; relaxed case solves coupled f and phi again. Field gradient, interaction, shell polarization, core monopole, shell monopole, exterior A and ledger recorded. Fixed-profile u=r delta_phi operator independently discretized at h=.05 with outer Dirichlet wall. All min eigenvalues positive (minimum .00380057); sufficient negative-potential bound ≤.36<1 including exterior. This does not qualify full coupled nonlinear stability; shell is prescribed/pinned, no shell mechanical mass or stress dynamics. Absolute A−3pi Q ≤4.94e-8, R48 controls included. Exterior has Z=1 and beta=0 so 1/r asymptotic remains (tiny core tail truncated).

| shell r=4–8 | fixed Q/Qbase | relaxed Q/Qbase |
|---|---:|---:|
| none |1|1|
| stiffness .1x |1|.909715|
| stiffness 10x |1|1.03428|
| beta=.1 |.436474|.453848|
| beta=1 |.0407991|.0438114|
| beta=5 |.00116053|.00125908|
| beta=-.01 |1.13522|1.12498|
| beta=-.03 |1.54339|1.48897|

Derivative-only shell changes phi_center for fixed source by +17.8% at Z=.1 or -5.95% at Z=10, no change in flux. Relaxed core changes Q via alpha rho phi and profile. Positive beta produces opposite induced shell source, reducing outward flux; beta=5 shell almost cancels source (~99.9% attenuation). Negative beta yields same-sign source and enhancement; verified stable fixed-profile bound for explored values, not arbitrary strong negative response. Exterior is Coulombic, unlike globally massive Yukawa field in 28C. The near-field screening need not destroy the far-field 1/r law. But a shell which breaks shift symmetry supplies new source content and does not make gravity universal. Pinned shell also breaks translation invariance, so translation inertia/virial values in relaxed shell rows are diagnostics only, not dressed-system inertia. Energy of actual shell material missing.

## Programme decision from D/E/F

D: a covariant vacuum can induce EH stiffness with suitable field content, but geometry and universal metric coupling are inputs. E: self-energy/dressing substantially changes the source but does not enforce an energy relation. F: shift-symmetric stiffness does not shield fixed flux, while explicit polarizing material can; no GR shielding result. The promising next move is a shared microscopic origin for metric response and matter stress, not tuning a secondary field's polarizability to match Newton. A useful further test is covariant matter + self-field coupling from a common metric action, comparing stress-integral source to complete dressed energy (including binding), with nonmetric controls; separately explore how effective metric symmetry can arise from collective vacuum excitations without being inserted. Neither branch yet derives G, spacetime from scratch, or real gravitation.
