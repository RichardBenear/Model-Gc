"""EX-GC37A stage S1: induced (Sakharov) spin-2 stiffness from one lattice scalar on Y1.

Lattice units h = 1. Scalar: E(k)^2 = sum_a K_a^2 + m^2, K_a = 2 sin(k_a/2) (Y1 / 34N dispersion).
External TT metric h_ij(x, tau) = eps cos(q z + Omega tau), q along z, two polarisations:
  cross: h_xy on plaquettes, linear vertex -h Gbar_x Gbar_y (34N placement: link gradients averaged to the plaquette,
         f_x = i K_x cos(k_y/2), f_y = i K_y cos(k_x/2));
  plus:  h_xx = -h_yy at sites, linear vertex (h/2)(-G_x^2 + G_y^2), f_a = i K_a.
Second-order (contact) terms come from sqrt(gamma) gamma^ab and 1/sqrt(gamma) expanded to O(h^2), each h^2 placed as an
average of squares (sum-preserving), so the contact term is independent of (Omega, q).
Euclidean effective action Gamma = (1/2) int h Pi h, so for h = eps cos: Gamma/(beta V) = eps^2 Pi / 4.
  Pi_contact = int (K_x^2 + K_y^2) / (2E)
  Pi_bubble(Omega, q) = -(1/2) int F(k, -k-q)^2 (E1+E2) / (2 E1 E2 ((E1+E2)^2 + Omega^2)),  E1 = E(k), E2 = E(k+q)
Readouts: A = d Pi/d Omega^2, B = d Pi/d q^2 at 0; c_g^2 = B/A; 1/G_ind = 16 pi A (EH: Pi = Q^2/(16 pi G)).
Checks: G3 exact diagonalisation of the static lattice Hamiltonian; G2 continuum Pauli-Villars control (c_g = 1).
Plan: spec/EX-GC37A_plan.md. Our construction, not LaViolette's Model G.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC37A"


def Kf(k):
    return 2 * np.sin(k / 2)


def grid(n, shift=0.5):
    k = 2 * np.pi * (np.arange(n) + shift) / n - np.pi
    return np.meshgrid(k, k, k, indexing="ij")


def F2(kx, ky, kz, q, pol, lattice=True):
    """F(k, k') with k' = -k - q (q along z). Returns F^2 (F is real)."""
    if lattice:
        Kx, Ky = Kf(kx), Kf(ky)
        Kxp, Kyp = Kf(-kx), Kf(-ky)
        cx, cy = np.cos(kx / 2), np.cos(ky / 2)
    else:
        Kx, Ky, Kxp, Kyp, cx, cy = kx, ky, -kx, -ky, 1.0, 1.0
    if pol == "cross":
        # f_x f_y' + f_x' f_y with f_x = i K_x cos(k_y/2), f_y = i K_y cos(k_x/2); k'_x = -k_x, k'_y = -k_y
        F = -(Kx * cy * Kyp * cx + Kxp * cy * Ky * cx)
    else:
        # f_x f_x' - f_y f_y' with f_a = i K_a
        F = -(Kx * Kxp) + (Ky * Kyp)
    return F ** 2


def energy(kx, ky, kz, m, lattice=True):
    if lattice:
        return np.sqrt(Kf(kx) ** 2 + Kf(ky) ** 2 + Kf(kz) ** 2 + m ** 2)
    return np.sqrt(kx ** 2 + ky ** 2 + kz ** 2 + m ** 2)


def bubble(kx, ky, kz, w, m, Om, q, pol, lattice=True):
    E1 = energy(kx, ky, kz, m, lattice)
    E2 = energy(kx, ky, kz + q, m, lattice)
    s = E1 + E2
    return -0.5 * np.sum(w * F2(kx, ky, kz, q, pol, lattice) * s / (2 * E1 * E2 * (s * s + Om * Om)))


def delta_bubble(kx, ky, kz, w, m, Om, q, pol, lattice=True):
    """Pi_b(Om, q) - Pi_b(0, 0), formed inside the integrand for accuracy."""
    E1 = energy(kx, ky, kz, m, lattice)
    E2 = energy(kx, ky, kz + q, m, lattice)
    s, s0 = E1 + E2, 2 * E1
    a = F2(kx, ky, kz, q, pol, lattice) * s / (2 * E1 * E2 * (s * s + Om * Om))
    b = F2(kx, ky, kz, 0.0, pol, lattice) * s0 / (2 * E1 * E1 * s0 * s0)
    return -0.5 * np.sum(w * (a - b))


def contact(kx, ky, kz, w, m, lattice=True):
    E = energy(kx, ky, kz, m, lattice)
    Kx2 = Kf(kx) ** 2 if lattice else kx ** 2
    Ky2 = Kf(ky) ** 2 if lattice else ky ** 2
    return np.sum(w * (Kx2 + Ky2) / (2 * E))


def coeffs(n, m, pol, steps=(0.04, 0.08)):
    kx, ky, kz = grid(n)
    w = 1.0 / n ** 3                       # int d^3k/(2pi)^3 over the zone = mean over the grid (h = 1)
    rows = {}
    for s in steps:
        rows[s] = (delta_bubble(kx, ky, kz, w, m, s, 0.0, pol), delta_bubble(kx, ky, kz, w, m, 0.0, s, pol))
    s1, s2 = steps
    # fit dPi = A s^2 + A4 s^4 from two steps
    def fit(i):
        y1, y2 = rows[s1][i], rows[s2][i]
        M = np.array([[s1 ** 2, s1 ** 4], [s2 ** 2, s2 ** 4]])
        return np.linalg.solve(M, [y1, y2])[0]
    A, B = fit(0), fit(1)
    pc = contact(kx, ky, kz, w, m)
    pb0 = bubble(kx, ky, kz, w, m, 0.0, 0.0, pol)
    return dict(n=n, m=m, pol=pol, A=float(A), B=float(B), cg2=float(B / A), inv_G=float(16 * np.pi * A),
                G_ind=float(1 / (16 * np.pi * A)), contact=float(pc), bubble0=float(pb0),
                mass_term=float(pc + pb0), mass_term_over_A=float((pc + pb0) / A))


def g3(n=24, m=0.3, q_index=2, pol="cross", eps=1e-3):
    """Exact static vacuum energy E0(eps) = (1/2) sum of normal-mode frequencies on an n^3 periodic lattice,
    h(z) = eps cos(q z), block-diagonal in (k_x, k_y). Pi(0, q) = 4 E2 / V with E0 = E0(0) + eps^2 E2."""
    q = 2 * np.pi * q_index / n
    k = 2 * np.pi * np.arange(n) / n
    z = np.arange(n)

    def E0(e):
        hz = e * np.cos(q * z)
        hz2 = hz ** 2
        tot = 0.0
        for kx in k:
            for ky in k:
                Kx2, Ky2 = Kf(kx) ** 2, Kf(ky) ** 2
                a = 1 + hz2 / 2                                  # 1/sqrt(gamma) kinetic (inverse mass)
                diag = (1 + hz2 / 2) * (Kx2 + Ky2) + (1 - hz2 / 2) * m ** 2
                if pol == "cross":
                    diag = diag - 2 * hz * Kf(kx) * Kf(ky) * np.cos(kx / 2) * np.cos(ky / 2)
                else:
                    diag = diag + (-hz) * Kx2 + hz * Ky2
                czz = 1 - (hz2 + np.roll(hz2, -1)) / 4          # link z+1/2 between z and z+1
                M = np.diag(diag + czz + np.roll(czz, 1))
                M -= np.diag(czz[:-1], 1) + np.diag(czz[:-1], -1)
                M[0, -1] -= czz[-1]
                M[-1, 0] -= czz[-1]
                sa = np.sqrt(a)
                ev = np.linalg.eigvalsh(sa[:, None] * M * sa[None, :])
                tot += 0.5 * np.sum(np.sqrt(np.maximum(ev, 0)))
        return tot

    e0, ep, em = E0(0.0), E0(eps), E0(-eps)
    E2 = (ep + em - 2 * e0) / (2 * eps ** 2)
    V = n ** 3
    exact = 4 * E2 / V
    # same quantity from the momentum formulas on the same discrete grid
    kk = np.meshgrid(k, k, k, indexing="ij")
    w = 1.0 / V
    formula = contact(*kk, w, m) + bubble(*kk, w, m, 0.0, q, pol)
    return dict(n=n, m=m, q=float(q), pol=pol, exact=float(exact), formula=float(formula),
                rel_diff=float(abs(exact - formula) / abs(formula)))


def g2(m=0.3, pv=(3.0, 6.0), pol="cross", kmax=60.0, nr=400, nang=48):
    """Continuum Pauli-Villars control: Pi_reg = sum_i c_i Pi(m_i), sum c = 0, sum c m^2 = 0.
    Lorentz invariance requires Pi(Om = s, q = 0) = Pi(Om = 0, q = s)."""
    M1, M2 = pv
    # 1 + c1 + c2 = 0 ; m^2 + c1 M1^2 + c2 M2^2 = 0
    c1, c2 = np.linalg.solve([[1, 1], [M1 ** 2, M2 ** 2]], [-1, -m ** 2])
    # cylindrical-symmetric integrand about z: use spherical coords (r, theta, phi); Gauss-Legendre in r on [0, kmax]
    xr, wr = np.polynomial.legendre.leggauss(nr)
    r = 0.5 * kmax * (xr + 1); wr = 0.5 * kmax * wr
    xt, wt = np.polynomial.legendre.leggauss(nang)
    ph = 2 * np.pi * (np.arange(2 * nang) + 0.5) / (2 * nang); wp = 2 * np.pi / (2 * nang)
    R, CT, PH = np.meshgrid(r, xt, ph, indexing="ij")
    ST = np.sqrt(1 - CT ** 2)
    kx, ky, kz = R * ST * np.cos(PH), R * ST * np.sin(PH), R * CT
    W = (wr[:, None, None] * wt[None, :, None] * wp) * R ** 2 / (2 * np.pi) ** 3
    out = {}
    for s in (0.05, 0.1, 0.2):
        tO = sum(c * delta_bubble(kx, ky, kz, W, mm, s, 0.0, pol, lattice=False) for c, mm in ((1, m), (c1, M1), (c2, M2)))
        tq = sum(c * delta_bubble(kx, ky, kz, W, mm, 0.0, s, pol, lattice=False) for c, mm in ((1, m), (c1, M1), (c2, M2)))
        out[str(s)] = dict(dPi_Omega=float(tO), dPi_q=float(tq), ratio=float(tq / tO))
    return dict(m=m, pv_masses=list(pv), pv_coeffs=[float(c1), float(c2)], kmax=kmax, steps=out)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = sys.argv[1]
    t0 = time.process_time()
    if part == "g3":
        res = [g3(pol=p, q_index=qi) for p in ("cross", "plus") for qi in (1, 3)]
    elif part == "g2":
        res = [g2(pol="cross", kmax=km, nr=nr) for km, nr in ((60.0, 300), (120.0, 600), (240.0, 1200))] + [g2(pol="plus", kmax=120.0, nr=600)]
    elif part == "s1":
        res = []
        for pol in ("cross", "plus"):
            for m in (0.1, 0.5, 1.0):
                for n in (64, 128, 192):
                    r = coeffs(n, m, pol)
                    res.append(r)
                    print(r, flush=True)
    res = dict(part=part, rows=res, cpu_seconds=time.process_time() - t0)
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))
