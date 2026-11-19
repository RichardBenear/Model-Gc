"""EX-GC07 (exploratory): 3D time simulation of the clocked Gc3 particle, spectral in angle.

Fields u = h0 + delta, delta expanded in complex spherical harmonics Y_lm (l <= LMAX, all m) times a radial
finite-volume grid (cell centres r_i, dr = 0.05, R = 30). Diffusion implicit per (species, l) (LapN / LapD radial
operators of ex008, plus -l(l+1)/r^2 and zero Dirichlet for l >= 1; G Dirichlet at l = 0); reactions explicit,
evaluated on a Gauss-Legendre (8) x Fourier (16) angular grid. IMEX Euler: fixed points are exact steady states.
Gc3: W_t = ... - kappa V W ; tau_V V_t = D_V lap V + W - V (G feeds X, k = 0.2). Order: X, Y, W, G, V.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import sph_harm_y

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc06_clock as c6  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC07")
LMAX, NT, NP = 4, 8, 16
LM = [(l, m) for l in range(LMAX + 1) for m in range(-l, l + 1)]


class Sim:
    def __init__(self, kappa, tauV, dv=0.1, R=30.0, dr=0.05, dt=0.004, seed=0, amp=1e-3):
        self.kappa, self.tauV, self.dt = kappa, tauV, dt
        M, u0 = steady(kappa, dv, R, dr)
        self.M, self.u0 = M, u0
        self.r, self.n = M.r, M.n
        x, wgl = np.polynomial.legendre.leggauss(NT)
        th = np.arccos(x)
        ph = 2 * np.pi * np.arange(NP) / NP
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        self.th, self.ph = TH.ravel(), PH.ravel()
        self.wq = (wgl[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]).ravel()
        self.Y = np.array([sph_harm_y(l, m, self.th, self.ph) for l, m in LM]).T     # (nq, nlm)
        self.Ya = (np.conj(self.Y) * self.wq[:, None]).T                              # analysis (nlm, nq)
        self.h0 = np.array([e.A, e.B / e.A, 0.0, e.GINF, 0.0])
        self.D = np.array([1.0, e.DY, e.DW, 1.0, dv])
        self.tau = np.array([1, 1, 1, 1, tauV], float)
        # delta initial: the steady radial particle in l = 0 (Y00 = 1/sqrt(4 pi))
        self.d = np.zeros((5, len(LM), self.n), complex)
        self.d[:, 0, :] = (u0 - self.h0[:, None]) * np.sqrt(4 * np.pi)
        rng = np.random.default_rng(seed)
        shape = u0[2] / u0[2].max()
        for i, (l, m) in enumerate(LM):
            if l <= 2:
                self.d[2, i, :] += amp * (rng.standard_normal() + 1j * rng.standard_normal()) * shape
        self._real_fix()
        cen = sp.diags(1.0 / self.r ** 2)
        self.lu = {}
        I = sp.identity(self.n, format="csc")
        for s in range(5):
            for l in range(LMAX + 1):
                if l == 0:
                    L = M.LapD if s == 3 else M.LapN
                else:
                    L = M.LapD - l * (l + 1) * cen
                self.lu[s, l] = spla.splu((I - dt / self.tau[s] * self.D[s] * L).tocsc())
        self.lidx = {l: [i for i, (ll, m) in enumerate(LM) if ll == l] for l in range(LMAX + 1)}
        self.shell = self.r ** 2 * dr
        X, Yc, Zc = (np.sin(self.th) * np.cos(self.ph), np.sin(self.th) * np.sin(self.ph), np.cos(self.th))
        self.dirs = np.array([X, Yc, Zc])                                             # (3, nq)
        self.t = 0.0

    def _real_fix(self):
        """Enforce the reality condition c_{l,-m} = (-1)^m conj(c_{l,m})."""
        for i, (l, m) in enumerate(LM):
            if m < 0:
                j = LM.index((l, -m))
                self.d[:, i, :] = (-1) ** m * np.conj(self.d[:, j, :])
            elif m == 0:
                self.d[:, i, :] = self.d[:, i, :].real

    def phys(self):
        return self.h0[:, None, None] + np.real(np.matmul(self.Y, self.d))

    def reaction(self, u):
        X, Yv, W, G, V = u
        s = c6.K * G * W
        return np.array([e.A - (1 + e.B) * X + X * X * Yv + s,
                         e.B * X - X * X * Yv,
                         -W + e.C * Yv * W * W - W ** 3 - self.kappa * V * W,
                         -s,
                         W - V])

    def step(self):
        u = self.phys()
        Rq = self.reaction(u) / self.tau[:, None, None]
        Rc = np.matmul(self.Ya, Rq.astype(complex))
        rhs = self.d + self.dt * Rc
        for s in range(5):
            for l, idx in self.lidx.items():
                b = rhs[s, idx, :].T
                k = b.shape[1]
                x = self.lu[s, l].solve(np.ascontiguousarray(np.hstack([b.real, b.imag])))
                self.d[s, idx, :] = (x[:, :k] + 1j * x[:, k:]).T
        self.t += self.dt

    def diag(self):
        u = self.phys()
        W = u[2]
        vol = self.wq[:, None] * self.shell[None, :]
        Wint = float(np.sum(vol * W))
        cen = [float(np.sum(vol * W * self.r[None, :] * self.dirs[a][:, None]) / Wint) for a in range(3)]
        Mc = float(np.sum(vol * c6.K * u[3] * W))
        amps = {f"l{l}": float(np.sqrt(np.sum(np.abs(self.d[2, self.lidx[l], :]) ** 2 * self.shell))) for l in
                range(1, LMAX + 1)}
        m1 = {f"m{m}": float(np.sqrt(np.sum(np.abs(self.d[2, LM.index((1, m)), :]) ** 2 * self.shell)))
              for m in (-1, 0, 1)}
        return dict(t=self.t, M=Mc, Wint=Wint, Wc=float(W[:, 0].mean()), cx=cen[0], cy=cen[1], cz=cen[2],
                    **amps, **m1, finite=bool(np.all(np.isfinite(u))))


def steady(kappa, dv, R, dr):
    z = np.load(os.path.join(e.ROOT, "results", "EX-GC06", f"particle_kappa{kappa:g}_dv{dv:g}.npz"))
    M = c6.Clock(kappa, dv, R=R, dr=dr)
    u = np.array([np.interp(M.r, z["r"], z["u"][s]) for s in range(5)])
    u, res, ok = M.newton(u)
    assert ok, res
    return M, u


def run(kappa, tauV, T, every=0.5, tag=None, **kw):
    os.makedirs(OUT, exist_ok=True)
    s = Sim(kappa, tauV, **kw)
    tag = tag or f"k{kappa:g}_tV{tauV:g}"
    rows, nstep, k_every = [], int(round(T / s.dt)), int(round(every / s.dt))
    t0 = time.time()
    for i in range(nstep + 1):
        if i % k_every == 0:
            rows.append(s.diag())
            if not rows[-1]["finite"]:
                break
        if i % (k_every * 50) == 0:
            r = rows[-1]
            print(tag, "t %.1f M %.5f Wc %.4f c (%.3f %.3f %.3f) l1 %.2e l2 %.2e  wall %.0f" % (
                r["t"], r["M"], r["Wc"], r["cx"], r["cy"], r["cz"], r["l1"], r["l2"], time.time() - t0), flush=True)
            json.dump(dict(kappa=kappa, tauV=tauV, dt=s.dt, rows=rows, wall=time.time() - t0),
                      open(os.path.join(OUT, f"ts_{tag}.json"), "w"))
        if i < nstep:
            s.step()
    json.dump(dict(kappa=kappa, tauV=tauV, dt=s.dt, rows=rows, wall=time.time() - t0),
              open(os.path.join(OUT, f"ts_{tag}.json"), "w"))
    np.savez_compressed(os.path.join(OUT, f"final_{tag}.npz"), d=s.d, r=s.r, t=s.t)
    print(tag, "done wall", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    kap, tV, T = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    extra = dict(a.split("=") for a in sys.argv[4:])
    run(kap, tV, T, tag=extra.pop("tag", None), **{k: float(v) if k != "seed" else int(v) for k, v in extra.items()})
