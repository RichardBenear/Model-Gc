"""EX-GC39B B2: can a lapse well hold its own matter? Complex scalar chi (mass mu) in -N^2 dt^2 + N^-2 dx^2,
mode chi = A u(r) e^{-i w t}:  lap u + (w^2 N^-4 - mu^2 N^-2) u = 0;  lap N = 4 pi G rho,  N(R) = 1,
rho = A^2 (w^2 N^-4 u^2 + u'^2 + mu^2 N^-2 u^2). Radial lattice h = 1 (the cutoff). Plan: spec/EX-GC39B_plan.md.
Our construction, not LaViolette's Model G.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC39B"


def lump(G, mu, M, h=1.0, R=4000.0, newtonian=False, tol=1e-11, maxit=4000, node=0, mix=0.3):
    n = int(R / h)
    r = h * np.arange(1, n + 1)
    N = np.ones(n)
    for it in range(maxit):
        if newtonian:     # Schroedinger-Newton: -(1/2mu) y'' + mu (N - 1) y = E y   (N - 1 = Phi)
            main = 1 / (mu * h * h) + mu * (N - 1); off = -0.5 / (mu * h * h) * np.ones(n - 1)
            ev, vec = eigh_tridiagonal(main, off, select="i", select_range=(node, node))
            y = vec[:, 0]; w2 = (mu + ev[0]) ** 2
            u = y / r
            dens = u * u
        else:
            Bh = N ** -2.0                       # B^{1/2}, B = N^-4
            main = (2 / h ** 2 + mu * mu * N ** -2) / Bh ** 2
            off = -1 / h ** 2 / (Bh[1:] * Bh[:-1])
            ev, vec = eigh_tridiagonal(main, off, select="i", select_range=(node, node))
            w2 = ev[0]
            y = vec[:, 0] / Bh
            u = y / r
            up = np.gradient(u, r)
            dens = w2 * N ** -4 * u * u + up * up + mu * mu * N ** -2 * u * u
        shell = 4 * np.pi * r * r * h
        A2 = M / np.sum(dens * shell)
        rho = A2 * dens
        m = np.cumsum(rho * shell)
        outer = np.cumsum((4 * np.pi * r * rho * h)[::-1])[::-1] - 4 * np.pi * r * rho * h
        Nn = 1 - G * (m / r + outer) + G * M / R
        if np.min(Nn) <= 0.05:
            return dict(ok=False, why="lapse collapse", it=it, N0=float(np.min(Nn)))
        diff = np.max(np.abs(Nn - N))
        N = (1 - mix) * N + mix * Nn
        if diff < tol:
            break
    frac = m / m[-1]
    R90 = float(np.interp(0.9, frac, r))
    w = float(np.sqrt(w2)) if w2 > 0 else float("nan")
    bound = w < mu if not newtonian else True
    return dict(ok=bool(diff < tol and bound), it=it, w=w, w_over_mu=w / mu, N0=float(N[0]), R90=R90,
                C=G * mu * mu * M * R90, GmuM=G * mu * M, residual=float(diff), r_edge_frac=float(R90 / R),
                nodes=int(np.sum(np.diff(np.sign(y[np.abs(y) > 1e-12 * np.max(np.abs(y))])) != 0)))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "edge":
    # descriptive: where the family ends (mu = 0.1, G = 1), slower mixing
    rows = []
    for gm in (0.5, 0.55, 0.6, 0.62, 0.64, 0.66, 0.7):
        r_ = lump(1.0, 0.1, gm / 0.1, mix=0.05, maxit=40000); r_["GmuM"] = gm; rows.append(r_)
        print(gm, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r_.items()}, flush=True)
    (OUT / "b2_edge.json").write_text(json.dumps(rows, indent=1, default=float) + "\n")
    sys.exit()

if __name__ == "__main__":
    t0 = time.process_time()
    res = dict(scan=[], checks={})
    for G, mu in ((1.0, 0.05), (1.0, 0.1), (1.0, 0.2), (21.0, 0.1)):
        for gm in (0.02, 0.04, 0.08, 0.16, 0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0):
            M = gm / (G * mu)
            row = dict(G=G, mu=mu, M=M, **lump(G, mu, M))
            row["GmuM"] = gm
            if row.get("ok") and row["R90"] > 1000:
                row["note"] = "R90 > 1000: box too small, excluded"
            res["scan"].append(row)
            print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    # separate path: continuum Schroedinger-Newton in scaled units G = mu = M = 1, fine grid
    sn = lump(1.0, 1.0, 1.0, h=0.005, R=60.0, newtonian=True)
    res["checks"]["SN_continuum_C"] = sn["C"]
    res["checks"]["SN_continuum_R90"] = sn["R90"]
    # G = 0: no bound state (w >= mu)
    z = lump(1e-12, 0.1, 1.0, R=2000.0, maxit=3)
    res["checks"]["G0_w_over_mu"] = z["w_over_mu"]
    # first radial excitation, one case
    res["checks"]["node1_G1_mu0.1_gm0.3"] = lump(1.0, 0.1, 3.0, node=1)
    res["cpu_seconds"] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "b2.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(res["checks"])
