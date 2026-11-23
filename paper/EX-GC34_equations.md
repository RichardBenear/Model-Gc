# The GC34 lattice equations

Working reference through EX-GC34M | 2026-11-23

This document collects the equations used in the GC34 studies and connects them to their implementations. It covers the static Regge benchmark, the J1 linear vacuum model, the prescribed moving-source benchmark, and the scalar matter-source test. GC34 is a reference model developed within the Model Gc project. Its equations are separate from the original reaction–diffusion equations for the consumed carrier \(G\).

## 1. Grid and variables

The spatial grid is cubic with spacing \(h\). Each cube is divided into six tetrahedra sharing a body diagonal. An edge of type \(e\) runs from vertex \(v\) to \(v+d_e\), where the seven forward offsets are

\[
(1,0,0),\ (0,1,0),\ (0,0,1),\
(1,1,0),\ (0,1,1),\ (1,0,1),\ (1,1,1).
\]

Write \(l_e^0=h|d_e|\), \(\hat d_e=d_e/|d_e|\), and \(l_e=l_e^0+\delta l_e\). The gravity variables at a vertex are the lapse perturbation \(n_v=N_v-1\), the three shift components \(\boldsymbol\beta_v\), and the seven outgoing edge-length perturbations. The lapse sets the local clock rate. The shift accounts for the motion of coordinate labels between time slices. Edge lengths determine the spatial geometry. The code uses \(c=1\) and \(\kappa=8\pi G\).

[Grid and edge definitions](../src/exgc34d.py)

## 2. Static gravity — GC34D

For an interior edge, the curvature measure is the deficit angle

\[
\epsilon_e=2\pi-\sum_{\tau\supset e}\theta_{\tau e},
\]

where \(\theta_{\tau e}\) is the interior dihedral angle at edge \(e\) in tetrahedron \(\tau\). Let \(\bar N_e=(N_v+N_{v+d_e})/2\). The static gravitational term and the energy-only source term are

\[
L_{\rm static}=\frac{1}{\kappa}\sum_e \bar N_e l_e\epsilon_e,
\qquad
L_{\rm source}=-\sum_v N_v M_v.
\]

Here \(M_v\) is energy per vertex, so \(M_v/h^3\) is the corresponding density. The static Fourier matrix \(C\) is the Hessian of \(\kappa L_{\rm static}\) at \(N=1,\ l_e=l_e^0\). The implementation obtains it by differentiating the tetrahedron contributions and assembling them across the lattice.

For the static, energy-only solution studied in GC34D,

\[
(\Delta_h n)_v=\frac{\kappa}{2h^3}(M_v-\langle M\rangle),
\qquad
\frac{\delta l_e}{l_e^0}=-\frac{n_v+n_{v+d_e}}{2}
\]

in the conformal representative, up to a lattice coordinate displacement. Here \(\Delta_h\) is the nearest-neighbor spatial Laplacian. The mean subtraction is the convention in the periodic box. In the far field of an isolated mass the continuum limit is \(n=-GM/r\). The accompanying spatial perturbation is the part of the benchmark that gives the GR value of weak light bending in the straight-ray calculation.

[Static kernel](../src/exgc34d.py) · [GC34D results](../results/EX-GC34D/EX-GC34D_note.md)

## 3. Time-dependent gravity — GC34E and J1

The edge strain-rate variable used in the kinetic term is

\[
\sigma_e(v)=l_e^0\left[
 \dot{\delta l}_e(v)-
 \hat d_e\cdot\big(\boldsymbol\beta(v+d_e)-\boldsymbol\beta(v)\big)
\right].
\]

The second term subtracts the change attributed to motion of the coordinate grid. For GC34E, the six edge rates in each tetrahedron determine a symmetric strain-rate tensor \(K_\tau\) via \(D_e^i D_e^j K_{\tau,ij}=\sigma_e\), with \(D_e=h d_e\). The quadratic kinetic term is

\[
L_{\rm kin}^{\rm E}
=\frac{1}{2\kappa}\sum_\tau V_\tau
\left[K_{\tau,ij}K_\tau^{ij}-(\operatorname{tr}K_\tau)^2\right].
\]

This is the construction whose vacuum spectrum contains the growing off-axis mode. J1 retains the edge strain-rate variables and uses a different quadratic form, determined by the re-slicing condition below.

### The equation implemented in the solver

Use Fourier dependence \(\exp(iq\cdot v-i\omega t)\), where \(q=h\mathbf k\) is dimensionless. For

\[
x=(n,\boldsymbol\beta,\delta l_1,\ldots,\delta l_7)^\mathsf T,
\]

the strain-rate map \(S\) is given by

\[
(Sx)_e=l_e^0\left[-i\omega\,\delta l_e-
 \big(e^{iq\cdot d_e}-1\big)\hat d_e\cdot\boldsymbol\beta\right].
\]

J1 uses a range-one Hermitian strain kernel \(\mathcal K(q)\). Its coefficients are obtained by a constrained fit to the GC34E kinetic kernel, with the continuum normalization \(\mu=-1\). The assembled equation is

\[
\boxed{
\mathcal M(q,\omega)x=f,\qquad
\mathcal M(q,\omega)=C(q)+S(q,\omega)^\dagger\mathcal K(q)S(q,\omega).
}
\]

The quadratic kinetic contribution to the action for the mode is \(\sigma^\dagger\mathcal K\sigma/(2\kappa)\). If \(\mathcal M=C_0+\omega B+\omega^2 A\), the equivalent time-domain equation at fixed \(q\) is

\[
C_0 x+iB\dot x-A\ddot x=f(t).
\]

In this polynomial \(C_0=\mathcal M(q,0)\) is the full zero-frequency matrix. It includes the shift block of the kinetic term. The static lapse-and-edge matrix in Section 2 is its restriction to \(\boldsymbol\beta=0\).

### Re-slicing symmetry and the vacuum test

J1 requires an exact null direction of \(\mathcal M\) for changes in time slicing:

\[
\mathcal M z_{\rm slice}=0,
\qquad
z_{\rm slice}=\big(-i\omega,-\mathbf g(q),0_7\big)^\mathsf T,
\qquad
g_a(q)=\frac{e^{iq_a}-1}{h}.
\]

This condition is imposed together with the three spatial displacement symmetries. The vacuum calculation sets \(f=0\) and finds physical wave modes after accounting for those freedoms. In synchronous gauge, \(n=\boldsymbol\beta=0\), the edge evolution is \(A_{ee}\ddot{\delta l}=C_{ee}\delta l\), with initial data satisfying the lapse and momentum constraints. J1 passed the registered vacuum evolution test in GC34K. That test used the sampled directions and spacings described in its plan.

[Original kinetic matrix](../src/exgc34e.py) · [J1 strain map and re-slicing identity](../src/exgc34i.py) · [J1 candidate and evolution](../src/exgc34j.py) · [GC34K results](../results/EX-GC34K/EX-GC34K_note.md)

## 4. Prescribed moving source — GC34L

In the moving-source benchmark the density and current are supplied to the solver. With a Fourier density amplitude \(\rho\) in energy-per-vertex units, velocity \(\mathbf v\), and \(\omega=\mathbf k\cdot\mathbf v\), GC34L used the link current

\[
j_a=\rho v_a e^{iq_a/2}\frac{q_a/2}{\sin(q_a/2)},
\qquad
f=\kappa(\rho,-\mathbf j,0_7)^\mathsf T.
\]

The ratio is defined by continuity at \(q_a=0\). This form obeys the re-slicing source condition \(z_{\rm slice}^\dagger f=0\). The spatial source vanishes in this benchmark at the order retained in velocity. GC34L also used a curvature readout invariant under the four linear gauge directions. The moving-source curvature error remained first order in spacing under that test.

[Moving-source implementation](../src/exgc34l.py) · [GC34L results](../results/EX-GC34L/EX-GC34L_note.md)

## 5. Dynamical scalar and its sources — GC34M

The GC34M scalar packet evolves on the flat background according to the discrete Klein–Gordon equation

\[
\boxed{
\ddot\phi_v=(\Delta_h\phi)_v-m^2\phi_v,\qquad
(\Delta_h\phi)_v=\frac{1}{h^2}
 \sum_{a=1}^{3}(\phi_{v+\hat a}+\phi_{v-\hat a}-2\phi_v).
}
\]

Here \(\hat a\) denotes one grid step on axis \(a\). The matter action uses tetrahedral gradient energy and volume assigned to vertices. For tetrahedron \(\tau\) with vertices \(0,1,2,3\), define

\[
(G_\tau)_{ab}=\tfrac12(l_{0a}^2+l_{0b}^2-l_{ab}^2),\qquad
V_\tau=\frac{\sqrt{\det G_\tau}}{6},\qquad
(\Delta\phi_\tau)_a=\phi_a-\phi_0 \quad(a,b=1,2,3).
\]

The gradient energy on that tetrahedron is \(\tfrac12V_\tau\,\Delta\phi_\tau^\mathsf T G_\tau^{-1}\Delta\phi_\tau\). Its volume contributes \(V_\tau/4\) to each vertex. With this vertex volume \(V_v\), the scalar time-derivative and mass terms in the Lagrangian are

\[
\frac12\sum_v V_v\big[(1-n_v)\dot\phi_v^2-(1+n_v)m^2\phi_v^2\big]
\]

at the linear order in geometry used for sourcing. The gradient energy enters with a minus sign, and its lapse coupling uses the average lapse perturbation on each link with the flat gradient weights. The shift couples to momentum using the forward axis links.

The matter variation gives the sources in the gravity equation:

\[
s_n=\frac{\partial L_{\rm matter}}{\partial n},\qquad
\mathbf s_\beta=\frac{\partial L_{\rm matter}}{\partial\boldsymbol\beta},\qquad
s_e=\frac{\partial L_{\rm matter}}{\partial\delta l_e},
\qquad
f=-\kappa(s_n,\mathbf s_\beta,s_e)^\mathsf T.
\]

In this convention \(-s_n\) is the energy per vertex. The flat-background shift variation is

\[
s_{\beta_a}(v)=-h^2\frac{\dot\phi_v+\dot\phi_{v+\hat a}}2
                   (\phi_{v+\hat a}-\phi_v).
\]

The two local balances tested by GC34M are

\[
\begin{aligned}
0&=-\dot s_n(v)+\sum_a\frac{s_{\beta_a}(v)-s_{\beta_a}(v-\hat a)}{h},\\
\mathbf 0&=-\dot{\mathbf s}_\beta(v)+
\sum_e\hat d_e\big[s_e(v-d_e)-s_e(v)\big].
\end{aligned}
\]

The first balance, associated with re-slicing, held to roundoff. The second, associated with spatial displacement, had a first-order local residual even though total momentum was conserved in the periodic box. GC34M evolved the scalar on a flat background and tested its candidate sources. A self-consistent scalar-plus-gravity evolution remains a separate experiment.

[GC34M scalar action and source code](../src/exgc34m.py) · [Registered plan](../spec/EX-GC34M_plan.md) · [GC34M results](../results/EX-GC34M/EX-GC34M_note.md)
