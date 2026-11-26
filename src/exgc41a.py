"""EX-GC41A: delay and bending from the arrival-time front on a random pair network (deterministic link times).

Plane-wave source (all sites x < -140), receiver strip x in [118, 122], paired flat/lens runs on each graph.
Continuum benchmark: a plane-wave ray family through n = 1/v(N) for the same source and receiver lines.
Plan: spec/EX-GC41_plan.md. Our construction, not LaViolette's Model G.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc37b as B37

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC41"
TAU = B37.TAU
HALF, RC, RS = 150.0, 1.5, 3.0
XSRC, XREC = -140.0, 120.0
LENSES = {"main": 2.0, "weak": 0.5}
R_CORE = {"main": 5.0, "weak": 0.0}
BS = (20, 40, 60)


def vrule(N, two):
    return B37.rate(N, N, True) if two else N


def dlnv_dN(N, two):
    if not two:
        return 1.0 / N
    e = np.exp(-N * TAU)
    v = 2 * N * (1 - e)
    return (2 * (1 - e) + 2 * N * TAU * e) / v


def lapse_xy(x, y, GM):
    return 1 - 2 * GM / np.sqrt(x * x + y * y + RS * RS)


# ---------------------------------------------------------------- continuum ray family
def ray_family(GM, two, y0=np.arange(8.0, 100.0, 0.05), dx=0.02, rcore=5.0):
    """Rays start horizontal at x = XSRC; y'' = (1+y'^2)(d_y ln n - y' d_x ln n), T = int n sqrt(1+y'^2) dx.
    Returns y at XREC, delay in equivalent length (T - T_flat) v(1), and the monotonicity flag."""
    v1 = vrule(1.0, two)
    def lnn_grad(x, y):
        rho2 = x * x + y * y + RS * RS
        N = 1 - 2 * GM / np.sqrt(rho2)
        g = -dlnv_dN(N, two)                    # d ln n / dN
        dN = 2 * GM / rho2 ** 1.5
        return g * dN * x, g * dN * y, N
    def f(x, z):
        y, p, T = z
        gx, gy, N = lnn_grad(x, y)
        n = 1.0 / vrule(N, two)
        return np.array([p, (1 + p * p) * (gy - p * gx), n * np.sqrt(1 + p * p)])
    z = np.array([y0.copy(), np.zeros_like(y0), np.zeros_like(y0)])
    nstep = int(round((XREC - XSRC) / dx))
    x = XSRC
    rmin = np.full(len(y0), np.inf)
    with np.errstate(all="ignore"):
        for _ in range(nstep):
            k1 = f(x, z); k2 = f(x + dx / 2, z + dx / 2 * k1); k3 = f(x + dx / 2, z + dx / 2 * k2); k4 = f(x + dx, z + dx * k3)
            z = z + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4); x += dx
            rmin = np.minimum(rmin, np.hypot(x, z[0]))
    yend, T = z[0], z[2]
    delay = (T - (XREC - XSRC) / v1) * v1
    # rays that enter the removed core (r < 5) or blow up are dropped, as the graph has no sites there
    valid = np.isfinite(yend) & np.isfinite(T) & (rmin >= rcore) & (np.abs(yend) < HALF)
    yv, dv = yend[valid], delay[valid]
    return yv, dv, bool(np.all(np.diff(yv) > 0))


def continuum(GM, two):
    yend, delay, mono = ray_family(GM, two)
    order = np.argsort(yend)
    ye, de = yend[order], delay[order]
    out = {"monotone": mono}
    for b in BS:
        m5 = (ye >= b - 5) & (ye <= b + 5); m10 = (ye >= b - 10) & (ye <= b + 10)
        yy = np.linspace(b - 5, b + 5, 201); D = float(np.mean(np.interp(yy, ye, de)))
        yy10 = np.linspace(b - 10, b + 10, 401); slope = np.polyfit(yy10, np.interp(yy10, ye, de), 1)[0]
        out[b] = dict(D=D, alpha=float(-slope), covered=bool(ye.min() < b - 10 and ye.max() > b + 10))
    return out


def born_alpha(GM, two, b, L=(XSRC, XREC)):
    """Thin-lens Born deflection along the straight line y = b: alpha = int d_y ln n dx (positive = towards the mass)."""
    x = np.linspace(L[0], L[1], 200001)
    rho2 = x * x + b * b + RS * RS
    N = 1 - 2 * GM / np.sqrt(rho2)
    dlnn_dy = -dlnv_dN(N, two) * 2 * GM / rho2 ** 1.5 * b
    return float(-np.trapezoid(dlnn_dy, x))


# ---------------------------------------------------------------- graph
def graph(seed):
    rng = np.random.default_rng([20261126, 41, seed])
    n = rng.poisson((2 * HALF) ** 2)
    X = (rng.random((n, 2)) * 2 - 1) * HALF
    pairs = cKDTree(X).query_pairs(RC, output_type="ndarray")
    return X, pairs


def adjacency(X, pairs, n, two, keep):
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    m = keep[i] & keep[j]; i, j = i[m], j[m]
    w = 1.0 / B37.rate(n[i], n[j], True) if two else 1.0 / n[i]
    return csr_matrix((w, (i, j)), shape=(len(X), len(X)))


def one_graph(seed):
    X, pairs = graph(seed)
    r = np.hypot(X[:, 0], X[:, 1])
    src = np.where(X[:, 0] < XSRC)[0]
    rec = np.where((X[:, 0] >= XREC - 2) & (X[:, 0] <= XREC + 2))[0]
    band = np.where((np.abs(X[:, 1]) < 80) & (X[:, 0] > -100) & (X[:, 0] < 140))[0]
    out = {}
    for two in (True, False):
        rule = "two" if two else "one"
        for lens, GM in LENSES.items():
            keep = r >= R_CORE[lens]
            A0 = adjacency(X, pairs, np.ones(len(X)), two, keep)
            A1 = adjacency(X, pairs, lapse_xy(X[:, 0], X[:, 1], GM), two, keep)
            T0 = dijkstra(A0, indices=src[keep[src]], min_only=True)
            T1 = dijkstra(A1, indices=src[keep[src]], min_only=True)
            bb = band[np.isfinite(T0[band])]
            vflat = 1.0 / np.polyfit(X[bb, 0], T0[bb], 1)[0]
            ok = rec[np.isfinite(T0[rec]) & np.isfinite(T1[rec])]
            y = X[ok, 1]; Dl = (T1[ok] - T0[ok]) * vflat; F0 = T0[ok] * vflat
            res = dict(vflat=float(vflat))
            for b in BS:
                m5 = np.abs(y - b) <= 5; m10 = np.abs(y - b) <= 10
                res[b] = dict(D=float(Dl[m5].mean()), alpha=float(-np.polyfit(y[m10], Dl[m10], 1)[0]),
                              null_slope=float(np.polyfit(y[m10], F0[m10], 1)[0]), n5=int(m5.sum()), n10=int(m10.sum()))
            out[f"{rule}_{lens}"] = res
    # P5 isotropy: point source at the centre, flat, two-site
    keep = np.ones(len(X), bool)
    A0 = adjacency(X, pairs, np.ones(len(X)), True, keep)
    c = int(np.argmin(np.hypot(X[:, 0], X[:, 1])))
    T = dijkstra(A0, indices=c)
    d = np.hypot(*(X - X[c]).T); m = (d > 95) & (d < 105) & np.isfinite(T)
    ang = np.arctan2(X[m, 1] - X[c, 1], X[m, 0] - X[c, 0]); k = ((ang + np.pi) / (2 * np.pi) * 16).astype(int) % 16
    v = d[m] / T[m]
    out["iso"] = [float(v[k == q].mean()) for q in range(16)] if m.sum() > 1000 else None
    return out


def boot(vals, nb=2000, seed=0):
    vals = np.asarray(vals); rng = np.random.default_rng(seed)
    bs = np.array([vals[rng.integers(0, len(vals), len(vals))].mean() for _ in range(nb)])
    return float(vals.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), float(bs.std())


def run_graphs(seeds):
    from multiprocessing import Pool
    import os
    with Pool(int(os.environ.get("MODELG_WORKERS", "4"))) as p:
        return p.map(one_graph, list(seeds))


def summarise(per, cont):
    S = {}
    for key in per[0]:
        if key == "iso":
            continue
        rule, lens = key.split("_")
        S[key] = {}
        for b in BS:
            c = cont[f"{rule}_{lens}"][b]
            D = boot([p[key][b]["D"] for p in per]); a = boot([p[key][b]["alpha"] for p in per])
            nl = boot([p[key][b]["null_slope"] for p in per])
            S[key][b] = dict(D=D, alpha=a, null=nl, D_c=c["D"], alpha_c=c["alpha"],
                             D_over_c=D[0] / c["D"], alpha_over_c=a[0] / c["alpha"],
                             null_halfwidth_over_alpha_c40=(nl[2] - nl[1]) / 2 / abs(cont[f"{rule}_{lens}"][40]["alpha"]))
    R = {}
    for lens in LENSES:
        for b in BS:
            t = [p[f"two_{lens}"][b]["alpha"] for p in per]; o = [p[f"one_{lens}"][b]["alpha"] for p in per]
            ratio = np.mean(t) / np.mean(o)
            rng = np.random.default_rng(1); bs = []
            for _ in range(2000):
                ii = rng.integers(0, len(t), len(t)); bs.append(np.mean(np.array(t)[ii]) / np.mean(np.array(o)[ii]))
            cr = cont[f"two_{lens}"][b]["alpha"] / cont[f"one_{lens}"][b]["alpha"]
            R[f"{lens}_{b}"] = dict(ratio=float(ratio), ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                                    continuum_ratio=float(cr), ratio_over_continuum=float(ratio / cr))
    iso = np.array([p["iso"] for p in per if p["iso"] is not None])
    secmean = iso.mean(0); se = iso.std(0, ddof=1).mean() / np.sqrt(len(iso)) / secmean.mean()
    I = dict(n=len(iso), spread=float((secmean.max() - secmean.min()) / secmean.mean()), se_rel=float(se),
             spread_over_se=float((secmean.max() - secmean.min()) / secmean.mean() / se))
    return dict(per_rule=S, ratios=R, isotropy=I)


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if part == "continuum":
        res = {}
        for two in (True, False):
            for lens, GM in LENSES.items():
                res[f"{'two' if two else 'one'}_{lens}"] = continuum(GM, two)
        res["H2_zero"] = {r: continuum(0.0, r == "two")[60] for r in ("two", "one")}
        res["H2_born_weak_b60"] = {r: dict(born=born_alpha(0.5, r == "two", 60), ray=res[f"{r}_weak"][60]["alpha"])
                                   for r in ("two", "one")}
        # descriptive, added after H2 was seen: Born limit at tiny GM, and step-size convergence
        tiny = {}
        for r in ("two", "one"):
            yb, db, _ = ray_family(0.01, r == "two", y0=np.arange(40, 80, 0.05)); o = np.argsort(yb); yy = np.linspace(50, 70, 401)
            tiny[r] = float(-np.polyfit(yy, np.interp(yy, yb[o], db[o]), 1)[0] / born_alpha(0.01, r == "two", 60))
        res["born_limit_GM0.01_b60_ray_over_born"] = tiny
        conv = {}
        for r in ("two", "one"):
            ye, de, _ = ray_family(2.0, r == "two", dx=0.01); o = np.argsort(ye); yy = np.linspace(30, 50, 401)
            conv[r] = float(-np.polyfit(yy, np.interp(yy, ye[o], de[o]), 1)[0])
        res["dx_half_alpha_main_40"] = conv
    elif part == "pilot":
        per = run_graphs(range(100, 105))
        cont = json.loads((OUT / "continuum.json").read_text())["result"]
        cont = {k: {int(b): v for b, v in d.items() if b not in ("monotone",)} for k, d in cont.items() if k.startswith(("two", "one"))}
        a = np.array([p["two_main"][40]["alpha"] for p in per])
        sd = float(a.std(ddof=1)); ac = cont["two_main"][40]["alpha"]
        Ng = int(min(60, max(10, np.ceil((1.96 * sd / (0.02 * abs(ac))) ** 2))))
        res = dict(per_graph=per, sd_alpha_two_main_40=sd, alpha_c=ac, N_g=Ng)
    elif part == "scored":
        Ng = json.loads((OUT / "pilot.json").read_text())["result"]["N_g"]
        per = run_graphs(range(200, 200 + Ng))
        cont = json.loads((OUT / "continuum.json").read_text())["result"]
        cont = {k: {int(b): v for b, v in d.items() if b not in ("monotone",)} for k, d in cont.items() if k.startswith(("two", "one"))}
        res = dict(N_g=Ng, per_graph=per, summary=summarise(per, cont))
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res["result"] if part != "scored" else res["result"]["summary"], indent=1, default=float)[:4000])
