"""EX-GC06 (exploratory): a Gc core with its own clock ("Gc3", our construction).

Gc1 (G feeds X, k = 0.2) plus a slow core inhibitor V of W, zero in the ether:
  W_t = D_W lap W - W + c Y W^2 - W^3 - kappa V W
  tau_V V_t = D_V lap V + W - V
The steady particle depends on kappa and D_V, not on tau_V. For angular index l, eigenvalues of J_l psi = lam T psi,
T = diag(1, 1, 1, 1, tau_V) (order X, Y, W, G, V).
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC06")
K = 0.2
SHIFTS = [0.0, 0.1j, 0.25j, 0.5j, 1j, 2j, 0.3, 0.3 + 0.3j, 1.0]   # positive shifts catch fast-growing modes


class Clock(e.Model):
    def __init__(self, kappa, dv=0.1, R=120.0, dr=0.05):
        super().__init__(R, dr, K, 1, 0)
        self.kappa, self.dv = kappa, dv

    def F(self, u):
        X, Y, W, G, V = u
        f4 = super().F(u[:4])
        f4[2] = f4[2] - self.kappa * V * W
        fv = self.dv * (self.LapN @ V) + W - V
        return np.vstack([f4, fv[None]])

    def J(self, u, l=0):
        X, Y, W, G, V = u
        n, d = self.n, sp.diags
        if l == 0:
            LN, LG = self.LapN, self.LapD
        else:
            cen = d(-l * (l + 1) / self.r ** 2)
            LN = LG = self.LapD + cen
        k, kap = self.k, self.kappa
        return sp.bmat([
            [LN + d(-(1 + e.B) + 2 * X * Y), d(X * X), d(k * G), d(k * W), None],
            [d(e.B - 2 * X * Y), e.DY * LN - d(X * X), None, None, None],
            [None, d(e.C * W * W), e.DW * LN + d(-1 + 2 * e.C * Y * W - 3 * W * W - kap * V), None, d(-kap * W)],
            [None, None, d(-k * G), LG + d(-k * W), None],
            [None, None, sp.identity(n), None, self.dv * LN - sp.identity(n)]]).tocsc()

    def newton(self, u, tol=1e-10, maxit=40):
        for it in range(maxit):
            Fv = self.F(u).ravel()
            res = float(np.max(np.abs(Fv)))
            if res < tol:
                return u, res, True
            u = u + spla.spsolve(self.J(u), -Fv).reshape(5, self.n)
            if not np.all(np.isfinite(u)) or np.any(u[[0, 1, 3]] <= 0):
                return u, np.inf, False
        res = float(np.max(np.abs(self.F(u))))
        return u, res, res < tol

    def consumption(self, u):
        return self.integral(self.k * u[3] * u[2])


def particle(kappa, dv=0.1, steps=8):
    z = np.load(os.path.join(e.ROOT, "results", "EX-GC01", "radial_ref.npz"))
    u = np.vstack([z["u"], z["u"][2][None]])
    for t in np.linspace(0, 1, steps + 1)[1:]:
        M = Clock(t * kappa, dv)
        u, res, ok = M.newton(u)
        if not ok or u[2, 0] < 0.5:
            return None, None
    return M, u


def spectrum(M, u, l, tauV, k=16, keep=-1.0):
    J = M.J(u, l)
    T = np.concatenate([np.ones(4 * M.n), np.full(M.n, tauV)])
    A = (sp.diags(1 / T) @ J).tocsc()
    core = np.tile(M.r < 10, 5) & ~((np.arange(5 * M.n) >= 3 * M.n) & (np.arange(5 * M.n) < 4 * M.n))
    found = []
    for sg in SHIFTS:
        try:
            vals, vecs = spla.eigs(A, k=k, sigma=sg + 1e-4, which="LM")
        except Exception as ex:  # noqa: BLE001
            print("eigs failed", l, sg, ex, flush=True)
            continue
        for lam, v in zip(vals, vecs.T):
            if lam.real < keep or lam.imag < -1e-9:
                continue
            if any(abs(lam - complex(f["re"], f["im"])) < 1e-7 * max(1, abs(lam)) for f in found):
                continue
            w = np.abs(v) ** 2
            found.append(dict(re=float(lam.real), im=float(lam.imag), core=float(w[core].sum() / w.sum())))
    found.sort(key=lambda f: -f["re"])
    return found


def leading(sp_, l):
    """Leading core mode, skipping the l = 1 translation mode and the slow G ladder."""
    for m in sp_:
        if m["core"] < 0.8:
            continue
        if l == 1 and abs(m["im"]) < 1e-6 and m["re"] > -2e-3:
            continue
        return m
    return None


def run(kappa, dv, tauVs, ls=(0, 1, 2, 3)):
    M, u = particle(kappa, dv)
    tag = f"kappa{kappa:g}_dv{dv:g}"
    if M is None:
        print(tag, "no particle", flush=True)
        json.dump(dict(kappa=kappa, dv=dv, particle=False), open(os.path.join(OUT, f"{tag}.json"), "w"))
        return
    np.savez_compressed(os.path.join(OUT, f"particle_{tag}.npz"), r=M.r, u=u)
    info = dict(kappa=kappa, dv=dv, particle=True, M=M.consumption(u), W0=float(u[2, 0]), V0=float(u[4, 0]),
                rows=[])
    print(tag, "M", info["M"], "W0", info["W0"], flush=True)
    for tV in tauVs:
        for l in ls:
            s = spectrum(M, u, l, tV)
            info["rows"].append(dict(tauV=tV, l=l, lead=leading(s, l), modes=s[:8]))
            m = info["rows"][-1]["lead"]
            print(tag, "tauV", tV, "l", l, "lead", None if m is None else f"{m['re']:+.5f}{m['im']:+.4f}i c{m['core']:.2f}",
                  flush=True)
        json.dump(info, open(os.path.join(OUT, f"{tag}.json"), "w"), indent=1)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    kap, dv = float(sys.argv[1]), float(sys.argv[2])
    run(kap, dv, [float(x) for x in sys.argv[3:]])
