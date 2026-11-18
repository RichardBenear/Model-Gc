"""EX-GC01 (exploratory): weighted action and reaction for unequal Model Gc particles.

Radial part: M and mu when one particle's core parameters differ (k, c, D_W).
Pair part: full-domain axisymmetric (rho, z) co-moving solve with two speeds v_a (particle at z = -d/2, type a) and
v_b (z = +d/2, type b); parameters are piecewise (a for z < 0, b for z > 0; W ~ 0 at the midpoint).
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
from ex009_pair2d import faces  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC01")
REF = dict(k=0.2, c=1.0, dw=0.1)          # MG195 P1 (G feeds X)


class RadialP(e.Model):
    """Radial model with per-instance c and D_W (G feeds X)."""
    def __init__(self, R, dr, k, c, dw):
        super().__init__(R, dr, k, 1, 0)
        self.c, self.dw = c, dw

    def F(self, u):
        f = super().F(u)
        X, Y, W, G = u
        f[2] = self.dw * (self.LapN @ W) - W + self.c * Y * W * W - W ** 3
        return f

    def J(self, u, l=0):
        X, Y, W, G = u
        Jm = super().J(u, l).tolil()
        n = self.n
        if l == 0:
            LW = self.LapN
        else:
            LW = self.LapD + sp.diags(-2.0 / self.r ** 2)
        blkWY = sp.diags(self.c * W * W)
        blkWW = self.dw * LW + sp.diags(-1 + 2 * self.c * Y * W - 3 * W * W)
        Jm[2 * n:3 * n, n:2 * n] = blkWY
        Jm[2 * n:3 * n, 2 * n:3 * n] = blkWW
        return Jm.tocsc()


def radial(k, c, dw, R=120.0, dr=0.05, path=None):
    """Continue from the MG195 P1 profile to (k, c, dw) in 8 steps."""
    prof = np.load(os.path.join(e.ROOT, "results", "MG195_inverse_square", "radial_P1_R120.npz"))
    u = prof["u"].copy()
    ok = True
    for t in np.linspace(0, 1, 9)[1:]:
        p = {q: REF[q] + t * (v - REF[q]) for q, v in dict(k=k, c=c, dw=dw).items()}
        M = RadialP(R, dr, p["k"], p["c"], p["dw"])
        u, res, ok = M.newton(u)
        if not ok or u[2, 0] < 0.5:
            return None
    out = dict(k=k, c=c, dw=dw, residual=res, M=M.integral(k * u[3] * u[2]), mu=M.drift(u)["mu"],
               W0=float(u[2, 0]), eig_l0=e.leading_eigs(M, u, 0)[0][0], eig_l1=e.leading_eigs(M, u, 1)[:2])
    if path:
        np.savez_compressed(path, r=M.r, u=u)
    return out


def scan():
    rows = []
    for k in (0.1, 0.15, 0.2, 0.3, 0.4):
        rows.append(radial(k, 1.0, 0.1))
    for c in (0.95, 0.975, 1.025, 1.05, 1.1):
        rows.append(radial(0.2, c, 0.1))
    for dw in (0.07, 0.085, 0.12, 0.15):
        rows.append(radial(0.2, 1.0, dw))
    return rows


if __name__ == "__main__" and sys.argv[1:2] == ["scan"]:
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rows = scan()
    for r in rows:
        print(r if r is None else {q: (round(v, 6) if isinstance(v, float) else v) for q, v in r.items()}, flush=True)
    json.dump(dict(rows=rows, wall=time.time() - t0), open(os.path.join(OUT, "radial_scan.json"), "w"), indent=1,
              default=float)


class AxiFull:
    """Full-domain axisymmetric grid, particles a at z = -d/2 and b at z = +d/2."""
    def __init__(self, d, pa, pb, h=0.1, rho_core=8.0, pad=8.0, far=150.0):
        self.d = d
        rf = faces(rho_core, h, far)
        zp = faces(d / 2 + pad, h, d / 2 + far)
        zf = np.concatenate([-zp[::-1], zp[1:]])
        rc, zc = 0.5 * (rf[1:] + rf[:-1]), 0.5 * (zf[1:] + zf[:-1])
        self.nr, self.nz = len(rc), len(zc)
        self.R, self.Z = np.meshgrid(rc, zc, indexing="ij")
        dr, dz = np.diff(rf), np.diff(zf)
        ring = np.pi * (rf[1:] ** 2 - rf[:-1] ** 2)
        self.V = (ring[:, None] * dz[None, :]).ravel()
        n = self.nr * self.nz
        self.n = n
        idx = np.arange(n).reshape(self.nr, self.nz)
        rows, cols, vals = [], [], []
        diag = np.zeros(n)
        for i in range(self.nr - 1):
            c = 2 * np.pi * rf[i + 1] * dz / (rc[i + 1] - rc[i])
            a, b = idx[i], idx[i + 1]
            rows += [a, b]; cols += [b, a]; vals += [c, c]
            diag[a] -= c; diag[b] -= c
        for j in range(self.nz - 1):
            c = ring / (zc[j + 1] - zc[j])
            a, b = idx[:, j], idx[:, j + 1]
            rows += [a, b]; cols += [b, a]; vals += [c, c]
            diag[a] -= c; diag[b] -= c
        diagD, bD = diag.copy(), np.zeros(n)
        cr = 2 * np.pi * rf[-1] * dz / (dr[-1] / 2)
        diagD[idx[-1]] -= cr; bD[idx[-1]] += cr
        for j, hz in ((0, dz[0]), (-1, dz[-1])):
            cz = ring / (hz / 2)
            diagD[idx[:, j]] -= cz; bD[idx[:, j]] += cz
        K0 = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
        iV = sp.diags(1 / self.V)
        self.LN = (iV @ (K0 + sp.diags(diag))).tocsr()
        self.LD = (iV @ (K0 + sp.diags(diagD))).tocsr()
        self.bD = bD / self.V
        rr, cc, vv = [], [], []
        for j in range(self.nz):
            a = idx[:, j]
            jl, jr = max(j - 1, 0), min(j + 1, self.nz - 1)
            hh = zc[jr] - zc[jl]
            rr += [a, a]; cc += [idx[:, jr], idx[:, jl]]; vv += [np.full(self.nr, 1 / hh), np.full(self.nr, -1 / hh)]
        self.Dz = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(n, n))
        Zf = self.Z.ravel()
        left = Zf < 0
        self.k = np.where(left, pa["k"], pb["k"])
        self.c = np.where(left, pa["c"], pb["c"])
        self.dw = np.where(left, pa["dw"], pb["dw"])
        self.wa = 0.5 * (1 - np.tanh(Zf / 2.0))
        self.wb = 1 - self.wa

    def F(self, u):
        X, Y, W, G = u
        s = self.k * G * W
        return np.array([self.LN @ X + e.A - (1 + e.B) * X + X * X * Y + s,
                         e.DY * (self.LN @ Y) + e.B * X - X * X * Y,
                         self.dw * (self.LN @ W) - W + self.c * Y * W * W - W ** 3,
                         self.LD @ G + self.bD * e.GINF - s])

    def J(self, u):
        X, Y, W, G = u
        d = sp.diags
        k = self.k
        return sp.bmat([
            [self.LN + d(-(1 + e.B) + 2 * X * Y), d(X * X), d(k * G), d(k * W)],
            [d(e.B - 2 * X * Y), e.DY * self.LN - d(X * X), None, None],
            [None, d(self.c * W * W), d(self.dw) @ self.LN + d(-1 + 2 * self.c * Y * W - 3 * W * W), None],
            [None, None, d(-k * G), self.LD + d(-k * W)]]).tocsc()

    def guess(self, ra, rb):
        """ra, rb: (r, u, M) radial solutions for particles a and b."""
        h0 = np.array([e.A, e.B / e.A, 0.0, e.GINF])
        out = np.tile(h0[:, None], (1, self.n))
        Rf, Zf = self.R.ravel(), self.Z.ravel()
        for (pr, pu, M), z0 in ((ra, -self.d / 2), (rb, self.d / 2)):
            r = np.sqrt(Rf ** 2 + (Zf - z0) ** 2)
            for s in range(3):
                out[s] += np.interp(r, pr, pu[s] - h0[s], right=0.0)
            near = np.interp(r, pr, pu[3] - h0[3])
            out[3] += np.where(r < 15, near, -M / (4 * np.pi * np.maximum(r, 1e-9)))
        return out

    def solve(self, u, tol=1e-9, maxit=25):
        n = self.n
        Dz4 = sp.block_diag([self.Dz] * 4).tocsr()
        WA, WB = np.tile(self.wa, 4), np.tile(self.wb, 4)
        ref = u.ravel().copy()
        dX = (Dz4 @ ref) * np.repeat([1.0, 0, 0, 0], n)
        pa, pb = dX * WA, dX * WB
        va = vb = 0.0
        for it in range(maxit):
            uf = u.ravel()
            du = Dz4 @ uf
            Fv = self.F(u).ravel() + va * WA * du + vb * WB * du
            g = np.array([pa @ (uf - ref), pb @ (uf - ref)])
            res = float(np.max(np.abs(Fv)))
            if res < tol and np.max(np.abs(g)) < 1e-9:
                return u, va, vb, res, True, it
            Jm = self.J(u) + sp.diags(va * WA + vb * WB) @ Dz4
            lu = spla.splu(Jm.tocsc())
            C = np.column_stack([WA * du, WB * du])
            y0 = lu.solve(-Fv)
            Y = np.column_stack([lu.solve(C[:, 0]), lu.solve(C[:, 1])])
            S = np.array([[pa @ Y[:, 0], pa @ Y[:, 1]], [pb @ Y[:, 0], pb @ Y[:, 1]]])
            dv = np.linalg.solve(S, np.array([pa @ y0, pb @ y0]) + g)
            sol = np.concatenate([y0 - Y @ dv, dv])
            u = u + sol[:-2].reshape(4, n)
            va, vb = va + sol[-2], vb + sol[-1]
        uf = u.ravel()
        du = Dz4 @ uf
        res = float(np.max(np.abs(self.F(u).ravel() + va * WA * du + vb * WB * du)))
        return u, va, vb, res, res < tol, maxit

    def integral_side(self, f, side):
        m = (self.Z.ravel() < 0) if side == "a" else (self.Z.ravel() > 0)
        return float(np.sum((self.V * f)[m]))


TYPES = {"ref": dict(k=0.2, c=1.0, dw=0.1), "dw0.2": dict(k=0.2, c=1.0, dw=0.2), "k0.4": dict(k=0.4, c=1.0, dw=0.1)}
PAIRS = [("ref", "ref", d) for d in (22.0, 30.0, 40.0)] + [("ref", "dw0.2", d) for d in (22.0, 30.0, 40.0)] + \
        [("ref", "k0.4", d) for d in (22.0, 30.0, 40.0)] + [("dw0.2", "ref", 30.0), ("k0.4", "ref", 30.0)]


def run_pairs():
    os.makedirs(os.path.join(OUT, "partial"), exist_ok=True)
    rad, prof = {}, {}
    for name, p in TYPES.items():
        path = os.path.join(OUT, f"radial_{name}.npz")
        rad[name] = radial(**p, path=path)
        z = np.load(path)
        prof[name] = (z["r"], z["u"], rad[name]["M"])
        print("radial", name, rad[name], flush=True)
    rows = []
    for ta, tb, d in PAIRS:
        t0 = time.time()
        a = AxiFull(d, TYPES[ta], TYPES[tb])
        u, va, vb, res, ok, it = a.solve(a.guess(prof[ta], prof[tb]))
        s = a.k * u[3] * u[2]
        row = dict(a=ta, b=tb, d=d, ok=bool(ok), residual=res, iters=it, va=va, vb=vb,
                   Ma_pair=a.integral_side(s, "a"), Mb_pair=a.integral_side(s, "b"),
                   pred_va=-rad[ta]["mu"] * rad[tb]["M"] / (4 * np.pi * d ** 2),    # +z = toward b
                   pred_vb=rad[tb]["mu"] * rad[ta]["M"] / (4 * np.pi * d ** 2), wall=time.time() - t0)
        row["Ma_va"] = rad[ta]["M"] * va
        row["Mb_vb"] = rad[tb]["M"] * vb
        row["momentum_ratio"] = -row["Ma_va"] / row["Mb_vb"]
        row["momentum_ratio_pred"] = rad[ta]["mu"] / rad[tb]["mu"]
        np.savez_compressed(os.path.join(OUT, "partial", f"pair_{ta}_{tb}_d{d:g}.npz"), u=u, va=va, vb=vb)
        print("pair", row, flush=True)
        rows.append(row)
        json.dump(dict(radial=rad, pairs=rows), open(os.path.join(OUT, "pairs.json"), "w"), indent=1, default=float)


if __name__ == "__main__" and sys.argv[1:2] == ["pairs"]:
    run_pairs()
