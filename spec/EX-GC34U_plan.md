# EX-GC34U — isolated reciprocal scalar gravity (EXPLORATORY)

2026-11-24 | ChatGPT | Richard authorized the next reciprocal test identified in GC34T. Cap 2 core-hours.

## Scope and equations

Use the spherically symmetric continuum Einstein–Klein–Gordon equations discretized in radius, as an isolated benchmark before a nonlinear Y1 implementation. This deliberately replaces the periodic standing packet with a spherical massive packet, preserving its positive total mass. Polar-areal metric: ds² = −α² dt² + a² dr² + r² dΩ². Define Φ = ∂r φ, Π = (a/α) ∂t φ, and a⁻² = 1 − 2λ M/r. Set κ=1, scalar mass=1, with λ multiplying the gravity coupling. Evolve φ and Π by the Klein–Gordon equation and solve the radial Einstein constraints at each RK4 stage:

- M' = r²[(1−2λM/r)(Π²+Φ²)+φ²]/4.
- (log α)' = a²[λM/r² + λ r p_r/2], with p_r = [(Π²+Φ²)/a²−φ²]/2.
- φdot = αΠ/a; Πdot = r⁻² ∂r(r² αΦ/a) − α a φ.
- Regular center M(0)=0; match the lapse at the outer radius to sqrt(1−2λM/R). The packet starts with φ=Π=0.4 exp(−r²/8) and remains negligible at the boundary by t=4. Flat control λ=0. Signed coupling values are used as a numerical derivative only.

Use R=24, h=.125,.0625,.03125, dt=h/4, times 1, 2, 4. R=48 at h=.0625 for the boundary check. Evolve λ=0, ±.01 and ±.005. The initial scalar data use Π held fixed across couplings. The radius is defined by sphere area and the time normalization by the exterior metric.

## Gates

- U0, isolated weak-field solution: positive total mass for positive coupling; max compactness 2λM/r < .05 at the largest positive coupling; relative outer mass drift <=.005 on the fine grid; boundary scalar amplitude <1e-8. Compare the profiles of M, φ and Π inside r<=8 between R=24 and 48 at h=.0625, λ=.01 and t=4, with relative L2 difference <=1e-5.
- U1, conservation before interpreting feedback: differentiate the radial mass constraint in the scalar evolution direction and solve it for Mdot by a separate radial quadrature. Compare against the Einstein flux Mdot = r² α Π Φ/(2 a³), an equation unused in evolution, in relative L2 at all scored times. Also compare centered coupling derivatives of the two sides using ±.005. At the fine spacing, raw and derivative residuals must be <=.03 and improve by at least 2x relative to the coarse grid at every scored time.
- U2, reciprocal signal: compare centered coupling derivatives of M from the self-gravitating runs and from applying the metric constraint to flat-control scalar snapshots. Their difference isolates the change in the gravity source due to updating the scalar. At t=4 its weighted L2 norm relative to the total zeroth-order mass times sqrt(8) on r<=8 must be >=1e-4. The signal profile at ±.005 and ±.01 must agree within 5% fine-grid relative L2, and the two finer spatial grids must agree within 5%. Record the signed change of enclosed mass coefficient at r=2. Compare the measured derivative of the full M profile against the derivative predicted by integrating the differentiated Hamiltonian constraint using the measured scalar response.

Overall PASS requires U0–U2. This tests a continuum nonlinear completion in spherical symmetry and provides a reference for reciprocal gravity. It does not certify a nonlinear Y1 action.
