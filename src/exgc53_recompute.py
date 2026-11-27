"""EX-GC53 checks C1 (CG Poisson vs DST for one shift component of the ring), C2 (53E dt/2, short), C3 (shift off reproduces 52).
EXPLORATORY."""
import json, sys, time, platform
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc53 as T, exgc52 as S


def c1(n=49):
    q, P = T.ring_sources(n, 2, a_r=8.0, s_r=2.5); J = T.angmom(P); kap = 1.0; B = T.shift_static(P, kap, J)
    a = 0; LP = T.link_pos(n); bnd = T.dipole_beta(J, LP[a], kap / (4 * np.pi))[a]
    m = n - 2; e = np.ones(m); D2 = sp.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]); Id = sp.identity(m)
    A = sp.kron(sp.kron(D2, Id), Id) + sp.kron(sp.kron(Id, D2), Id) + sp.kron(sp.kron(Id, Id), D2)
    u = bnd.copy(); u[1:-1, 1:-1, 1:-1] = 0; bc = S.lap_aniso(u, 0.0)
    b = (T.LAM * kap * P[a])[1:-1, 1:-1, 1:-1] - bc; x, flag = cg(-A.tocsr(), -b.ravel(), rtol=1e-13, maxiter=10000)
    ref = B[a][1:-1, 1:-1, 1:-1].ravel(); return dict(rel=float(np.abs(x - ref).max() / np.abs(ref).max()), flag=int(flag))


def c2():
    r1 = T.run_pair(0.2, T=8.0, dt=0.1); r2 = T.run_pair(0.2, T=8.0, dt=0.05)
    t1, P1 = np.array(r1["t"]), np.array(r1["P1"]); P2 = np.array(r2["P1"])[::2]
    sel = t1 >= 2 - 1e-9; return dict(max_rel=float(np.max(np.abs(P2[sel] - P1[sel]) / np.abs(P1[sel]))), kappa=[r1["kappa"], r2["kappa"]])


def c3():
    r = T.run_pair(0.2, T=2.0, shift_on=False); D = json.loads((S.OUT / "D.json").read_text())["result"]["0.2"]
    P52 = np.array(D["P1"][:len(r["P1"])]); return dict(max_abs=float(np.max(np.abs(np.array(r["P1"]) - P52))), scale=float(np.max(np.abs(P52))), kappa=[r["kappa"], D["kappa"]])


if __name__ == "__main__":
    out = {}; t0 = time.time()
    for k in sys.argv[1:]:
        out[k] = {"1": c1, "2": c2, "3": c3}[k](); print("C" + k, out[k], flush=True)
    (T.OUT / ("checks" + "".join(sys.argv[1:]) + ".json")).write_text(json.dumps(dict(study="EX-GC53", platform=platform.platform(), wall_s=time.time() - t0, result=out), indent=1, default=float))
