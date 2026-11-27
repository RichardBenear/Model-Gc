"""EX-GC52 checks. C1: moving solver at V = 0 vs exgc51.static. C2: linear anisotropic solve vs scipy CG on the same discrete
problem. C3: 52D packet-1 force at V = 0.2 with dt/2 (compared with D.json by exgc52_summary). EXPLORATORY."""
import json, sys, time, platform
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc52 as S, exgc51 as G


def c1(n=65):
    c = n // 2; psi, _, _ = S.packet_k(n, S.M, 0.0, (c, c, c), S.SIG, (0, 0, 0))
    N1, i1 = S.moving([psi], [S.M], 0.0, 1.0, n, (c, c, c), 0.0); N2, i2 = G.static([psi], [S.M], 0.0, 1.0, n, (c, c, c))
    return dict(maxdiff=float(np.abs(N1 - N2).max()), it_moving=i1, it_static=i2)


def c2(n=41, V=0.3, kap=1e-6):
    c = n // 2; k0 = S.k_for_v(V, S.M, 0.0); psi, _, _ = S.packet_k(n, S.M, 0.0, (c, c, c), 4.0, (k0, 0, 0))
    q = G.charge(psi, np.ones((n, n, n)), S.M, 0.0)
    N, info = S.moving([psi], [S.M], 0.0, kap, n, (c, c, c), V, qfix=q, gauss=False, tol=1e-16, maxit=60)
    m = n - 2; e = np.ones(m); D2 = sp.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]); Id = sp.identity(m)
    A = (1 - V * V) * sp.kron(sp.kron(D2, Id), Id) + sp.kron(sp.kron(Id, D2), Id) + sp.kron(sp.kron(Id, Id), D2)
    b = kap * q[1:-1, 1:-1, 1:-1].ravel(); x, flag = cg(-A.tocsr(), -b, rtol=1e-13, maxiter=5000)
    df = G.fN(N)[1:-1, 1:-1, 1:-1].ravel() + 0.5             # f - f(1): linear response (nonlinear O(kappa^2) ~ 1e-12 here)
    return dict(rel=float(np.abs(df - x).max() / np.abs(x).max()), cg_flag=int(flag), info=info)


def c3():
    return S.part_D(Vs=(0.2,), dt=0.05)


if __name__ == "__main__":
    which = sys.argv[1:] or ["12", "3"]; out = {}; t0 = time.time()
    if "12" in which:
        out["C1"] = c1(); print("C1", out["C1"], flush=True); out["C2"] = c2(); print("C2", out["C2"], flush=True)
        (S.OUT / "checks12.json").write_text(json.dumps(dict(study="EX-GC52", part="checks12", platform=platform.platform(), wall_s=time.time() - t0, result=out), indent=1, default=float))
    if "3" in which:
        out = dict(C3=c3()); (S.OUT / "checks3.json").write_text(json.dumps(dict(study="EX-GC52", part="checks3", platform=platform.platform(), wall_s=time.time() - t0, result=out), indent=1, default=float))
