"""EX-GC26A (exploratory): a conservative 3D topological soliton (the plain Skyrme model, no pion mass).

Field phi = (phi0, phi1, phi2, phi3), |phi| = 1 (an SU(2) element), vacuum (1, 0, 0, 0) at the box faces.
Skyrme units, c = 1. With D_ij = d_i phi . d_j phi the static energy density is
    e = tr D + (1/2) [ (tr D)^2 - D_ij D_ij ],          E / (12 pi^2) >= |B|  (Faddeev-Bogomolny bound)
and the Lorentz-invariant Lagrangian density is  K - e  with  K = phi_t.phi_t (1 + tr D) - sum_i (phi_t . d_i phi)^2.
Space: average of forward- and backward-difference versions (no zero-energy grid mode); the discrete forces are exact derivatives of the discrete energy and
Lagrangian. Relaxation: arrested Newton flow. Dynamics: the Euler-Lagrange equations of the discrete Lagrangian,
a 4x4 solve per site, with the |phi| = 1 constraint by a Lagrange multiplier; RK4 in time.
Separate path: radial hedgehog energy by shooting (radial()), sharing no code with the 3D solver.
Our construction (choice of model, discretisation and tests); the Skyrme model is Skyrme (1961).
"""
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.ndimage import gaussian_filter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC26A")
NORM = 12 * np.pi ** 2


# ---------------------------------------------------------------- separate path: radial hedgehog
def radial(R=40.0):
    """Shoot f(0) = pi, f'(0) = -a to f(inf) = 0; energy E/(12 pi^2) = (1/3pi) int r^2f'^2 + 2 s^2 (1 + f'^2) + s^4/r^2."""
    def rhs(r, y):
        f, g = y
        s, s2 = np.sin(f), np.sin(2 * f)
        return [g, (-2 * r * g - s2 * g * g + s2 + s * s * s2 / (r * r)) / (r * r + 2 * s * s)]

    def shoot(a):
        r0 = 1e-4
        ev = lambda r, y: y[0]
        ev.terminal = True
        ev2 = lambda r, y: y[1]
        ev2.terminal = True
        sol = solve_ivp(rhs, (r0, R), [np.pi - a * r0, -a], events=(ev, ev2), rtol=1e-11, atol=1e-13)
        if sol.t_events[0].size:
            return 1, sol           # crossed zero: a too large
        return -1, sol              # turned back up: a too small

    lo, hi = 1.0, 3.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if shoot(mid)[0] > 0:
            hi = mid
        else:
            lo = mid
    a = (lo + hi) / 2
    _, sol = shoot(a)
    # use the trajectory up to where it is still accurate (before it departs), then a k/r^2 tail
    r, f, g = sol.t, sol.y[0], sol.y[1]
    cut = np.argmax(r > 6.0)
    rr, ff, gg = r[:cut], f[:cut], g[:cut]
    k = ff[-1] * rr[-1] ** 2
    dens = lambda r_, f_, g_: r_ ** 2 * g_ ** 2 + 2 * np.sin(f_) ** 2 * (1 + g_ ** 2) + np.sin(f_) ** 4 / r_ ** 2
    core = np.trapezoid(dens(rr, ff, gg), rr)
    tail = quad(lambda x: dens(x, k / x ** 2, -2 * k / x ** 3), rr[-1], np.inf, limit=200)[0]
    E2 = np.trapezoid(rr ** 2 * gg ** 2 + 2 * np.sin(ff) ** 2, rr) + quad(
        lambda x: x * x * (2 * k / x ** 3) ** 2 + 2 * np.sin(k / x ** 2) ** 2, rr[-1], np.inf)[0]
    E = (core + tail) / (3 * np.pi)
    return dict(a=a, k=k, E=E, E2_frac=E2 / (3 * np.pi) / E, profile=dict(r=rr[::20].tolist(), f=ff[::20].tolist()))


# ---------------------------------------------------------------- 3D lattice
# Derivatives: the energy (and Lagrangian) is the average of an all-forward-difference and an all-backward-difference
# version. Each one-sided version sees a site-to-site (checkerboard) oscillation as a large gradient, so there is no
# zero-energy grid mode (plain central differences have one, and on them the knot escapes through it); the average is
# second-order accurate.
class Grid:
    def __init__(self, h, Lb):
        self.h = h
        n = int(round(2 * Lb / h)) + 1
        self.n = n
        ax = (np.arange(n) - (n - 1) / 2) * h
        self.ax = ax
        self.X = np.meshgrid(ax, ax, ax, indexing="ij")
        self.r = np.sqrt(sum(x * x for x in self.X))
        m = np.zeros((n, n, n), bool)
        m[1:-1, 1:-1, 1:-1] = True
        self.inner = m                    # sites that may move; the outer face is fixed vacuum
        mp = np.zeros((n, n, n), bool)
        mp[:-1, :-1, :-1] = True
        mm = np.zeros((n, n, n), bool)
        mm[1:, 1:, 1:] = True
        self.mask = {1: mp, -1: mm}

    def d(self, f, i, s):
        if s == 1:
            return (np.roll(f, -1, axis=i + 1) - f) / self.h * self.mask[1]
        return (f - np.roll(f, 1, axis=i + 1)) / self.h * self.mask[-1]

    def dT(self, P, i, s):
        """Adjoint of d(., i, s): sum_x P(x).d phi(x) -> gradient w.r.t. phi(y)."""
        P = P * self.mask[s]
        if s == 1:
            return (np.roll(P, 1, axis=i + 1) - P) / self.h
        return (P - np.roll(P, -1, axis=i + 1)) / self.h

    def strain(self, phi, s):
        dp = [self.d(phi, i, s) for i in range(3)]
        D = np.array([[(dp[i] * dp[j]).sum(0) for j in range(3)] for i in range(3)])
        return dp, D

    def energy_parts(self, phi):
        out = np.zeros(2)
        for s in (1, -1):
            dp, D = self.strain(phi, s)
            tr = D[0, 0] + D[1, 1] + D[2, 2]
            out += 0.5 * np.array([tr.sum(), 0.5 * (tr * tr - (D * D).sum((0, 1))).sum()])
        out *= self.h ** 3 / NORM
        return float(out[0]), float(out[1])

    def energy(self, phi):
        return sum(self.energy_parts(phi))

    def grad(self, phi):
        """d(sum e)/d phi per site (energy = h^3 sum e)."""
        g = 0.0
        for s in (1, -1):
            dp, D = self.strain(phi, s)
            tr = D[0, 0] + D[1, 1] + D[2, 2]
            for i in range(3):
                P = 2 * dp[i] * (1 + tr) - 2 * sum(D[i, j] * dp[j] for j in range(3))
                g = g + 0.5 * self.dT(P, i, s)
        return g

    def degree(self, phi):
        tot = 0.0
        for s in (1, -1):
            dp = [self.d(phi, i, s) for i in range(3)]
            Mx = np.stack([phi, dp[0], dp[1], dp[2]], 0)
            tot += 0.5 * np.linalg.det(np.moveaxis(Mx, (0, 1), (-1, -2))).sum()
        return float(tot * self.h ** 3 / (2 * np.pi ** 2))

    def vacuum(self):
        p = np.zeros((4,) + (self.n,) * 3)
        p[0] = 1
        return p

    def hedgehog(self, prof, centre=(0, 0, 0), scale=(1, 1, 1)):
        X = [(self.X[i] - centre[i]) / scale[i] for i in range(3)]
        r = np.sqrt(sum(x * x for x in X))
        f = prof(r)
        p = np.empty((4,) + r.shape)
        p[0] = np.cos(f)
        for i in range(3):
            p[i + 1] = np.sin(f) * X[i] / np.maximum(r, 1e-12)
        return self.fix(p)

    def fix(self, p):
        p = p / np.sqrt((p * p).sum(0))
        p[:, ~self.inner] = self.vacuum()[:, ~self.inner]
        return p


def tangent(phi, v):
    return v - phi * (phi * v).sum(0)


def relax(G, phi, dt_fac=0.15, maxit=20000, tol=1e-9, check=50):
    """Arrested Newton flow: phi'' = -grad (tangent), velocity zeroed whenever the energy rises."""
    dt = dt_fac * G.h
    v = np.zeros_like(phi)
    E = G.energy(phi)
    hist = [E]
    it = 0
    for it in range(1, maxit + 1):
        a = tangent(phi, -G.grad(phi)) * G.inner
        v = tangent(phi, v + dt * a)
        phi = G.fix(phi + dt * v)
        if it % 5 == 0:
            En = G.energy(phi)
            if En > E:
                v[:] = 0
            E = En
        if it % check == 0:
            hist.append(E)
            if len(hist) > 4 and abs(hist[-1] - hist[-3]) < tol * abs(hist[-1]):
                break
    return phi, dict(iters=it, E=G.energy(phi), hist=hist[::10])


# ---------------------------------------------------------------- Lorentz-invariant dynamics
def accel(G, phi, v):
    vv = (v * v).sum(0)
    GQ = 0.0
    GP = 0.0
    Mdv = 0.0
    M = 0.0
    for s in (1, -1):
        dp, D = G.strain(phi, s)
        dv = [G.d(v, i, s) for i in range(3)]
        tr = D[0, 0] + D[1, 1] + D[2, 2]
        vd = [(v * dp[i]).sum(0) for i in range(3)]
        for i in range(3):
            GQ = GQ + 0.5 * G.dT(2 * vv * dp[i] - 2 * vd[i] * v, i, s)
            GP = GP + 0.5 * G.dT(2 * dp[i] * (1 + tr) - 2 * sum(D[i, j] * dp[j] for j in range(3)), i, s)
        sdd = sum((dp[i] * dv[i]).sum(0) for i in range(3))
        Mdv = Mdv + 0.5 * (2 * sdd * v - sum(dv[i] * vd[i] + dp[i] * (dv[i] * v).sum(0) for i in range(3)))
        M = M + 0.5 * ((1 + tr)[None, None] * np.eye(4)[:, :, None, None, None]
                       - sum(dp[i][:, None] * dp[i][None, :] for i in range(3)))
    R = 0.5 * (GQ - GP) - Mdv
    Mt = np.moveaxis(M, (0, 1), (-2, -1))
    rhs = np.stack([np.moveaxis(R, 0, -1), np.moveaxis(phi, 0, -1)], -1)
    sol = np.linalg.solve(Mt, rhs)
    a = np.moveaxis(sol[..., 0], -1, 0)
    b = np.moveaxis(sol[..., 1], -1, 0)
    mu = (-vv - (phi * a).sum(0)) / (phi * b).sum(0)
    return (a + mu * b) * G.inner


def hamiltonian(G, phi, v):
    H = 0.0
    P = np.zeros(3)
    for s in (1, -1):
        dp, D = G.strain(phi, s)
        tr = D[0, 0] + D[1, 1] + D[2, 2]
        K = (v * v).sum(0) * (1 + tr) - sum((v * dp[i]).sum(0) ** 2 for i in range(3))
        e = tr + 0.5 * (tr * tr - (D * D).sum((0, 1)))
        H += 0.5 * float(((K + e) * G.mask[s]).sum())
        pi = 2 * ((1 + tr) * v - sum(dp[i] * (v * dp[i]).sum(0) for i in range(3)))
        P += 0.5 * np.array([-float((pi * dp[i]).sum()) for i in range(3)])
    return H * G.h ** 3 / NORM, (P * G.h ** 3 / NORM).tolist()


def evolve(G, phi, v, T, dt, every):
    rows = []
    n = int(round(T / dt))
    ke = int(round(every / dt))
    t0 = time.time()
    for i in range(n + 1):
        if i % ke == 0:
            H, P = hamiltonian(G, phi, v)
            rows.append(dict(t=round(i * dt, 4), H=H, P=P, B=G.degree(phi), wall=round(time.time() - t0, 1)))
            print("  t", rows[-1]["t"], "H", round(H, 7), "P", np.round(P, 7).tolist(), "B", round(rows[-1]["B"], 5),
                  flush=True)
        if i == n:
            break
        k1v = accel(G, phi, v)
        k1p = v
        k2v = accel(G, G.fix(phi + 0.5 * dt * k1p), v + 0.5 * dt * k1v)
        k2p = v + 0.5 * dt * k1v
        k3v = accel(G, G.fix(phi + 0.5 * dt * k2p), v + 0.5 * dt * k2v)
        k3p = v + 0.5 * dt * k2v
        k4v = accel(G, G.fix(phi + dt * k3p), v + dt * k3v)
        k4p = v + dt * k3v
        phi = G.fix(phi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))
        v = tangent(phi, v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)) * G.inner
    return phi, v, rows


PROF = lambda r: np.pi * np.exp(-r / 1.2)


def case_static(h, Lb, kind):
    G = Grid(h, Lb)
    if kind == "hedgehog":
        phi = G.hedgehog(PROF)
    elif kind == "lopsided":
        rng = np.random.default_rng(26)
        phi = G.hedgehog(PROF, centre=(0.3, -0.2, 0.1), scale=(1.4, 0.8, 1.0))
        noise = np.array([gaussian_filter(rng.normal(size=phi.shape[1:]), 2.0) for _ in range(4)])
        phi = G.fix(phi + 0.15 * noise / noise.std())
    elif kind == "control_deg0":
        phi = G.hedgehog(lambda r: 2.5 * r * np.exp(-r))     # f(0) = 0: no winding, same size bump
    B0, E0 = G.degree(phi), G.energy(phi)
    phi, info = relax(G, phi)
    E2, E4 = G.energy_parts(phi)
    return G, phi, dict(case=kind, h=h, Lb=Lb, n=G.n, B_start=B0, E_start=E0, B=G.degree(phi), E=info["E"],
                        E2=E2, E4=E4, iters=info["iters"])


def main(which):
    os.makedirs(OUT, exist_ok=True)
    t0 = time.process_time()
    if which == "radial":
        rec = radial()
        print("radial a", rec["a"], "E/12pi^2", rec["E"], "E2 frac", rec["E2_frac"])
    elif which.startswith("static"):
        _, h, Lb, kind = which.split(":")
        G, phi, rec = case_static(float(h), float(Lb), kind)
        print("static", rec)
        if kind == "hedgehog" and abs(float(h) - 0.2) < 1e-9:
            np.save(os.path.join(OUT, f"relaxed_h{h}_L{Lb}.npy"), phi)
    elif which.startswith("kick"):
        _, h, Lb, T = which.split(":")
        h, Lb, T = float(h), float(Lb), float(T)
        G = Grid(h, Lb)
        path = os.path.join(OUT, f"relaxed_h{h:g}_L{Lb:g}.npy")
        phi = np.load(path)
        E0 = G.energy(phi)
        rng = np.random.default_rng(2026)
        env = np.exp(-(G.r / 2.0) ** 2)
        dv = np.array([gaussian_filter(rng.normal(size=phi.shape[1:]), 3.0) for _ in range(4)]) * env
        v = tangent(phi, dv) * G.inner
        H, _ = hamiltonian(G, phi, v)
        v *= np.sqrt(0.05 * E0 / (H - E0))                    # kinetic energy = 5% of E0
        phi1, v1, rows = evolve(G, phi, v, T, 0.1 * h, 0.5)
        Eend_static = G.energy(phi1)
        phi2, info = relax(G, phi1)
        rec = dict(E0=E0, rows=rows, E_static_end=Eend_static, recovered_E=info["E"], recovered_B=G.degree(phi2),
                   recover_iters=info["iters"], centre_shift=None)
        print("kick: recovered E", info["E"], "B", rec["recovered_B"])
    rec.update(which=which, cpu_s=time.process_time() - t0, platform=sys.version.split()[0] + " numpy " + np.__version__)
    tag = which.replace(":", "_")
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
