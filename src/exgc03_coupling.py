"""EX-GC03 (exploratory): which coupling of consumed G makes free fall universal across particle types?

Couplings: consumed G becomes X (GX), Y (GY) or the core species W (GW). Types: k, c, D_W variants.
Per type: gravitational charge p = <psi+, P>, drag weights c_s = <psi+_s, U'_s>, Jordan-chain mass m(tau).
Only ratios are used (psi+ normalisation cancels): mu = -p / (c . tau), acceleration a = p / m.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc01_unequal as g1  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC03")
COUPLINGS = {"GX": (1, 0, 0), "GY": (0, 1, 0), "GW": (0, 0, 1)}
TYPES = {"ref": dict(k=0.2, c=1.0, dw=0.1), "k0.1": dict(k=0.1, c=1.0, dw=0.1), "k0.4": dict(k=0.4, c=1.0, dw=0.1),
         "c0.95": dict(k=0.2, c=0.95, dw=0.1), "c1.1": dict(k=0.2, c=1.1, dw=0.1),
         "dw0.05": dict(k=0.2, c=1.0, dw=0.05), "dw0.2": dict(k=0.2, c=1.0, dw=0.2)}


class RadialC(g1.RadialP):
    def __init__(self, R, dr, k, c, dw, cpl):
        super().__init__(R, dr, k, c, dw)
        self.cx, self.cy, self.cw = cpl

    def F(self, u):
        f = super().F(u)              # includes cx, cy terms via e.Model
        f[2] += self.cw * self.k * u[3] * u[2]
        return f

    def J(self, u, l=0):
        Jm = super().J(u, l).tolil()
        n = self.n
        X, Y, W, G = u
        Jm[2 * n:3 * n, 2 * n:3 * n] = Jm[2 * n:3 * n, 2 * n:3 * n] + sp.diags(self.cw * self.k * G)
        Jm[2 * n:3 * n, 3 * n:4 * n] = sp.diags(self.cw * self.k * W)
        return Jm.tocsc()


def forcing(M, u):
    X, Y, W, G = u
    P = np.zeros((4, M.n))
    P[0], P[1], P[2], P[3] = M.cx * M.k * W * M.r, M.cy * M.k * W * M.r, M.cw * M.k * W * M.r, -M.k * W * M.r
    return P.ravel()


def nulls(M, u):
    J1 = M.J(u, 1)
    Up = np.array([np.gradient(u[s], M.r) for s in range(4)]).ravel()
    lu, luT = spla.splu(J1), spla.splu(J1.T.tocsc())
    rv, lv = Up.copy(), np.random.default_rng(0).standard_normal(J1.shape[0])
    for _ in range(10):
        rv = lu.solve(rv); rv /= np.linalg.norm(rv)
        lv = luT.solve(lv); lv /= np.linalg.norm(lv)
    rv *= (Up @ Up) / (rv @ Up)
    return J1, lv, rv


def mass_raw(J1, lv, rv, Tvec):
    gam = float(lv @ (Tvec * rv))
    alpha = gam / float(lv @ rv)
    B = sp.bmat([[J1, sp.csr_matrix(rv[:, None])], [sp.csr_matrix(lv[None, :]), None]]).tocsc()
    w = spla.spsolve(B, np.concatenate([-(Tvec * rv - alpha * rv), [0.0]]))[:-1]
    return -float(lv @ (Tvec * w)), gam


def particle(cpl, p, R=120.0, dr=0.05):
    """Continue from the reference GX particle: first coupling (cx,cy,cw), then type parameters."""
    z = np.load(os.path.join(e.ROOT, "results", "EX-GC01", "radial_ref.npz"))
    u = z["u"].copy()
    ref = g1.REF
    for t in np.linspace(0, 1, 9)[1:]:
        c3 = tuple((1 - t) * a + t * b for a, b in zip((1, 0, 0), cpl))
        M = RadialC(R, dr, ref["k"], ref["c"], ref["dw"], c3)
        u, res, ok = M.newton(u)
        if not ok or u[2, 0] < 0.5:
            return None, None
    for t in np.linspace(0, 1, 9)[1:]:
        q = {s: ref[s] + t * (p[s] - ref[s]) for s in ("k", "c", "dw")}
        M = RadialC(R, dr, q["k"], q["c"], q["dw"], cpl)
        u, res, ok = M.newton(u)
        if not ok or u[2, 0] < 0.5:
            return None, None
    return M, u


def analyse_type(M, u):
    J1, lv, rv = nulls(M, u)
    n = M.n
    c = np.array([lv[s * n:(s + 1) * n] @ rv[s * n:(s + 1) * n] for s in range(4)])
    p = float(lv @ forcing(M, u))
    Mc = M.integral(M.k * u[3] * u[2])
    eig1 = e.leading_eigs(M, u, 1)[:2]
    eig0 = e.leading_eigs(M, u, 0)[0]
    return dict(c=c, p=p, M=Mc, J1=J1, lv=lv, rv=rv, eig_l0=eig0, eig_l1=eig1)


def run():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    out = {}
    store = {}
    for cname, cpl in COUPLINGS.items():
        out[cname] = {}
        for tname, p in TYPES.items():
            M, u = particle(cpl, p)
            if M is None:
                out[cname][tname] = None
                print(cname, tname, "no particle", flush=True)
                continue
            a = analyse_type(M, u)
            store[(cname, tname)] = (M, u, a)
            sgn = np.sign(a["c"].sum())
            out[cname][tname] = dict(c=(a["c"] * sgn).tolist(), p=a["p"] * sgn, M=a["M"], mu1=-a["p"] / a["c"].sum(),
                                     eig_l0=a["eig_l0"], eig_l1=a["eig_l1"])
            print(cname, tname, "M %.4f mu1 %.5f c %s" % (a["M"], out[cname][tname]["mu1"],
                  np.round(np.array(out[cname][tname]["c"]) / abs(a["c"].sum()), 4)), round(time.time() - t0), flush=True)
    # analyses per coupling
    summary = {}
    for cname in COUPLINGS:
        names = [t for t in TYPES if out[cname][t] is not None]
        C = np.array([out[cname][t]["c"] for t in names])          # rows: types, cols: X Y W G
        P = np.array([out[cname][t]["p"] for t in names])
        mu1 = -P / C.sum(1)
        spreadA = float(np.std(mu1) / abs(np.mean(mu1)))
        # (B) best response times: tau in [0.1, 30]^4, ether-stable tau_Y/tau_X < 8.64; minimise mu spread
        def spread(lt):
            tau = np.exp(lt)
            g = C @ tau
            if np.any(g <= 0):
                return 10.0
            mu = -P / g
            return float(np.std(mu) / abs(np.mean(mu)))
        cons = [{"type": "ineq", "fun": lambda lt: np.log(8.0) - (lt[1] - lt[0])}]
        best = None
        for start in ([0, 0, 0, 0], [0, 1, 2, 0], [0, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 3], [-1, 0, 2, 1]):
            r = minimize(spread, np.array(start, float), method="SLSQP", bounds=[(np.log(0.1), np.log(30))] * 4,
                         constraints=cons)
            if best is None or r.fun < best.fun:
                best = r
        tauB = np.exp(best.x)
        # (C) low-drag regime: tau minimising max_i gamma_i(tau)/gamma_i(1) at fixed |log tau| budget (tau_X = 1)
        def maxdrag(lt):
            tau = np.exp(np.concatenate([[0.0], lt]))
            return float(np.max((C @ tau) / C.sum(1)))
        cons2 = [{"type": "ineq", "fun": lambda lt: np.log(8.0) - lt[0]}]
        bestC = None
        for start in ([0, 0, 0], [1, 2, 0], [2, 3, 0], [0, 3, 0], [0, 0, 3], [2, 3, 3]):
            r = minimize(maxdrag, np.array(start, float), method="SLSQP", bounds=[(np.log(0.1), np.log(60))] * 3,
                         constraints=cons2)
            if bestC is None or r.fun < bestC.fun:
                bestC = r
        tauC = np.exp(np.concatenate([[0.0], bestC.x]))
        acc = {}
        for t in names:
            M, u, a = store[(cname, t)]
            m, g = mass_raw(a["J1"], a["lv"], a["rv"], np.repeat(tauC, M.n))
            acc[t] = dict(m_over_p=m / a["p"], gamma_over_gamma1=g / a["c"].sum(), mu=-a["p"] / g)
        mop = np.array([acc[t]["m_over_p"] for t in names])
        summary[cname] = dict(types=names, mu_tau1=dict(zip(names, mu1.tolist())), spreadA=spreadA,
                              best_tau=tauB.tolist(), spreadB=float(best.fun),
                              lowdrag_tau=tauC.tolist(), lowdrag_maxdrag=float(bestC.fun), lowdrag=acc,
                              spreadC_m_over_p=float(np.std(mop) / abs(np.mean(mop))))
        print("SUMMARY", cname, "A spread %.4f | B spread %.4f at tau %s | C maxdrag %.3f at tau %s, m/p spread %.4f" % (
            spreadA, best.fun, np.round(tauB, 2), bestC.fun, np.round(tauC, 2), summary[cname]["spreadC_m_over_p"]), flush=True)
    json.dump(dict(per_type=out, summary=summary, wall=time.time() - t0), open(os.path.join(OUT, "results.json"), "w"),
              indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))


if __name__ == "__main__":
    run()
