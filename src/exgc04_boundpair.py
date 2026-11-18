"""EX-GC04 (exploratory): does a bound pair of identical Gc particles fall like a single particle?

AxiN: full-domain axisymmetric (rho, z) grid with N particles on the axis at given centres, one co-moving speed per
particle (bordered Newton by block elimination), and an external G gradient: outer boundary G = G_inf + E z.
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

OUT = os.path.join(e.ROOT, "results", "EX-GC04")
K = 0.2


def stretch(lo, hi, h, far, grow=1.1):
    mid = np.arange(lo, hi + 1e-9, h)
    right = [mid[-1]]
    w = h
    while right[-1] < hi + far:
        w *= grow
        right.append(right[-1] + w)
    left = [mid[0]]
    w = h
    while left[-1] > lo - far:
        w *= grow
        left.append(left[-1] - w)
    return np.concatenate([np.array(left[::-1]), mid[1:-1], np.array(right)])


class AxiN:
    def __init__(self, centres, E=0.0, h=0.1, rho_core=8.0, pad=8.0, far=150.0):
        self.centres = np.array(centres, float)
        self.E = E
        rf = faces(rho_core, h, far)
        zf = stretch(self.centres.min() - pad, self.centres.max() + pad, h, far)
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
        diagD = diag.copy()
        gb = np.zeros(n)                       # sum over boundary faces of coeff * G_boundary
        cr = 2 * np.pi * rf[-1] * dz / (dr[-1] / 2)
        diagD[idx[-1]] -= cr
        gb[idx[-1]] += cr * (e.GINF + E * zc)
        for j, zface, hz in ((0, zf[0], dz[0]), (-1, zf[-1], dz[-1])):
            cz = ring / (hz / 2)
            diagD[idx[:, j]] -= cz
            gb[idx[:, j]] += cz * (e.GINF + E * zface)
        K0 = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
        iV = sp.diags(1 / self.V)
        self.LN = (iV @ (K0 + sp.diags(diag))).tocsr()
        self.LD = (iV @ (K0 + sp.diags(diagD))).tocsr()
        self.gb = gb / self.V
        rr, cc, vv = [], [], []
        for j in range(self.nz):
            a = idx[:, j]
            jl, jr = max(j - 1, 0), min(j + 1, self.nz - 1)
            hh = zc[jr] - zc[jl]
            rr += [a, a]; cc += [idx[:, jr], idx[:, jl]]; vv += [np.full(self.nr, 1 / hh), np.full(self.nr, -1 / hh)]
        self.Dz = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(n, n))
        Zf = self.Z.ravel()
        # Voronoi switches along z (tanh, width 2) for each particle
        cs = np.sort(self.centres)
        mids = 0.5 * (cs[1:] + cs[:-1])
        edges = np.concatenate([[-np.inf], mids, [np.inf]])
        self.w = []
        for i in range(len(cs)):
            lo = 1.0 if np.isinf(edges[i]) else 0.5 * (1 + np.tanh((Zf - edges[i]) / 2.0))
            hi = 1.0 if np.isinf(edges[i + 1]) else 0.5 * (1 - np.tanh((Zf - edges[i + 1]) / 2.0))
            self.w.append(np.broadcast_to(lo * hi, Zf.shape).astype(float))

    def F(self, u):
        X, Y, W, G = u
        s = K * G * W
        return np.array([self.LN @ X + e.A - (1 + e.B) * X + X * X * Y + s,
                         e.DY * (self.LN @ Y) + e.B * X - X * X * Y,
                         e.DW * (self.LN @ W) - W + e.C * Y * W * W - W ** 3,
                         self.LD @ G + self.gb - s])

    def J(self, u):
        X, Y, W, G = u
        d = sp.diags
        return sp.bmat([
            [self.LN + d(-(1 + e.B) + 2 * X * Y), d(X * X), d(K * G), d(K * W)],
            [d(e.B - 2 * X * Y), e.DY * self.LN - d(X * X), None, None],
            [None, d(e.C * W * W), e.DW * self.LN + d(-1 + 2 * e.C * Y * W - 3 * W * W), None],
            [None, None, d(-K * G), self.LD + d(-K * W)]]).tocsc()

    def guess(self, pr, pu, M):
        h0 = np.array([e.A, e.B / e.A, 0.0, e.GINF])
        out = np.tile(h0[:, None], (1, self.n))
        Rf, Zf = self.R.ravel(), self.Z.ravel()
        out[3] += self.E * Zf
        for z0 in self.centres:
            r = np.sqrt(Rf ** 2 + (Zf - z0) ** 2)
            for s in range(3):
                out[s] += np.interp(r, pr, pu[s] - h0[s], right=0.0)
            near = np.interp(r, pr, pu[3] - h0[3])
            out[3] += np.where(r < 15, near, -M / (4 * np.pi * np.maximum(r, 1e-9)))
        return out

    def solve(self, u, ref=None, tol=1e-10, maxit=25):
        n, N = self.n, len(self.centres)
        Dz4 = sp.block_diag([self.Dz] * 4).tocsr()
        Ws = [np.tile(w, 4) for w in self.w]
        ref = u.ravel().copy() if ref is None else ref.ravel().copy()
        dX = (Dz4 @ ref) * np.repeat([1.0, 0, 0, 0], n)
        phis = [dX * Wi for Wi in Ws]
        v = np.zeros(N)
        for it in range(maxit):
            uf = u.ravel()
            du = Dz4 @ uf
            Fv = self.F(u).ravel() + sum(v[i] * Ws[i] * du for i in range(N))
            g = np.array([ph @ (uf - ref) for ph in phis])
            res = float(np.max(np.abs(Fv)))
            if res < tol and np.max(np.abs(g)) < 1e-9:
                return u, v, res, True, it
            Jm = self.J(u) + sp.diags(sum(v[i] * Ws[i] for i in range(N))) @ Dz4
            lu = spla.splu(Jm.tocsc())
            Cm = np.column_stack([Ws[i] * du for i in range(N)])
            y0 = lu.solve(-Fv)
            Y = np.column_stack([lu.solve(Cm[:, i]) for i in range(N)])
            S = np.array([[phis[a] @ Y[:, b] for b in range(N)] for a in range(N)])
            dv = np.linalg.solve(S, np.array([ph @ y0 for ph in phis]) + g)
            u = u + (y0 - Y @ dv).reshape(4, n)
            v = v + dv
        uf = u.ravel()
        du = Dz4 @ uf
        res = float(np.max(np.abs(self.F(u).ravel() + sum(v[i] * Ws[i] * du for i in range(N)))))
        return u, v, res, res < tol, maxit

    def consumption(self, u, i):
        return float(np.sum(self.V * K * u[3] * u[2] * self.w[i]))


# ---------------------------------------------------------------- drivers (EX-GC04 runs)
PROF = os.path.join(e.ROOT, "results", "MG195_inverse_square", "radial_P1_R120.npz")
M_P1 = 5.317804244620182          # radial (R = 120) consumption of the MG195 P1 particle (k = 0.2, G feeds X)
EG = 1e-4                         # external gradient used for the mobility (central difference in E)


def _prof():
    p = np.load(PROF)
    return p["r"], p["u"]


def run_vd(ds):
    """Pair speed v(d) at E = 0 (v of the right particle; > 0 = repulsion)."""
    pr, pu = _prof()
    rows = []
    for d in ds:
        a = AxiN([-d / 2, d / 2])
        u, v, res, ok, it = a.solve(a.guess(pr, pu, M_P1), maxit=30)
        rows.append(dict(d=d, v_right=float(v[1]), v_left=float(v[0]), ok=bool(ok), res=res,
                         Mc=[a.consumption(u, 0), a.consumption(u, 1)]))
        print("d", d, "v", v[1], ok, flush=True)
    json.dump(rows, open(os.path.join(OUT, f"vd_axin_{ds[0]:g}.json"), "w"), indent=1)


def run_single(h=0.1, ginf=None):
    """Single particle in G = G_inf + E z, E = +-EG. Optional reduced background G_inf."""
    pr, pu = _prof()
    pu = pu.copy()
    if ginf is not None:
        e.GINF = ginf
        pu[3] = pu[3] - 1 + ginf
    out = dict(Ginf=e.GINF, h=h, runs={})
    for E in (EG, -EG):
        a = AxiN([0.0], E=E, h=h)
        u, v, res, ok, it = a.solve(a.guess(pr, pu, M_P1))
        out["runs"][str(E)] = dict(v=float(v[0]), res=res, ok=bool(ok), it=it, Mc=a.consumption(u, 0), n=a.n)
    out["mu"] = (out["runs"][str(EG)]["v"] - out["runs"][str(-EG)]["v"]) / (2 * EG)
    name = f"single_gradient_h{h:g}.json" if ginf is None else f"single_lowGinf_h{h:g}.json"
    json.dump(out, open(os.path.join(OUT, name), "w"), indent=1)
    print(name, out["mu"], flush=True)


def run_pair(d, h=0.1):
    """Pair at +-d/2: E = 0 (continuation start), then E = +-EG. mu_pair from the mean speed."""
    pr, pu = _prof()
    out = dict(d=d, runs={})
    prev = None
    for E in (0.0, EG, -EG):
        a = AxiN([-d / 2, d / 2], E=E, h=h)
        u0 = a.guess(pr, pu, M_P1) if prev is None else prev + np.array([0, 0, 0, 1.0])[:, None] * E * a.Z.ravel()[None, :]
        u, v, res, ok, it = a.solve(u0)
        if E == 0.0:
            prev = u
        out["runs"][str(E)] = dict(v=[float(x) for x in v], res=res, ok=bool(ok), it=it,
                                   Mc=[a.consumption(u, 0), a.consumption(u, 1)], n=a.n)
    vp, vm = out["runs"][str(EG)]["v"], out["runs"][str(-EG)]["v"]
    out["mu_pair"] = (np.mean(vp) - np.mean(vm)) / (2 * EG)
    out["tidal"] = ((vp[1] - vp[0]) - (vm[1] - vm[0])) / (2 * EG)
    out["M_pair"] = sum(out["runs"]["0.0"]["Mc"])
    json.dump(out, open(os.path.join(OUT, f"pair_gradient_d{d:g}_h{h:g}.json"), "w"), indent=1)
    print("pair", d, h, out["mu_pair"], out["M_pair"], flush=True)


if __name__ == "__main__":
    # python src/exgc04_boundpair.py vd 6 6.5 ... | single [h] | lowg | pair d [h]
    os.makedirs(OUT, exist_ok=True)
    cmd, args = sys.argv[1], [float(x) for x in sys.argv[2:]]
    t0 = time.time()
    if cmd == "vd":
        run_vd(args)
    elif cmd == "single":
        run_single(*args)
    elif cmd == "lowg":
        run_single(0.1, ginf=1 - M_P1 / (4 * np.pi * 7.145))
    elif cmd == "pair":
        run_pair(*args)
    print("wall s", round(time.time() - t0), flush=True)
