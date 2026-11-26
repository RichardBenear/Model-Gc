"""EX-GC40B: 37B event counting on an irregular pair network (2D Poisson sprinkling, links within r_c).

Sites fire at rate N(x). Two-site rule: a signal crosses a link at a coincidence (37B, tau = 0.01). One-site control:
a signal crosses when the sending site fires. Every link costs the same in flat space, whatever its length: the
geometry is counts. Plan: spec/EX-GC40_plan.md. Our construction; not LaViolette's Model G.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc37b as B37

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC40B"
TAU = B37.TAU
R_DISK, RC, GM, RS = 150.0, 1.5, 2.0, 3.0
R_CORE = 5.0          # deviation: the registered GM, r_s give N <= 0 for r < 2.65; sites with r < 5 are removed
NREAL = 20


def graph(seed, R=R_DISK, rc=RC):
    rng = np.random.default_rng([20261126, seed])
    n = rng.poisson(np.pi * R * R)
    r = R * np.sqrt(rng.random(n)); th = 2 * np.pi * rng.random(n)
    X = np.stack([r * np.cos(th), r * np.sin(th)], 1)
    pairs = cKDTree(X).query_pairs(rc, output_type="ndarray")
    return X, pairs


def lapse(X, mass=True):
    if not mass:
        return np.ones(len(X))
    r = np.hypot(X[:, 0], X[:, 1])
    return 1 - 2 * GM / np.sqrt(r * r + RS * RS)


def det_weights(n, i, j, two):
    """Mean crossing time i -> j (37B rate sum)."""
    return 1.0 / B37.rate(n[i], n[j], True) if two else 1.0 / n[i]


@njit(cache=True)
def hop_samples(ni, nj, two, tau, seed):
    """One sample of the 37B crossing time per directed link (same rule as exgc37b.walk for a single hop)."""
    np.random.seed(seed)
    out = np.empty(len(ni))
    for e in range(len(ni)):
        nx, ny = ni[e], nj[e]
        if not two:
            out[e] = np.random.exponential(1.0 / nx)
            continue
        t = 0.0
        tx = -np.random.exponential(1.0 / nx)
        ty = -np.random.exponential(1.0 / ny)
        while True:
            t += np.random.exponential(1.0 / (nx + ny))
            if np.random.random() < nx / (nx + ny):
                if t - ty <= tau:
                    break
                tx = t
            else:
                if t - tx <= tau:
                    break
                ty = t
        out[e] = t
    return out


def build(X, pairs, n, two, stochastic, seed, keep=None):
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    if keep is not None:
        m = keep[i] & keep[j]; i, j = i[m], j[m]
    if stochastic:
        w = hop_samples(n[i].astype(float), n[j].astype(float), two, TAU, seed)
    else:
        w = det_weights(n, i, j, two)
    return csr_matrix((w, (i, j)), shape=(len(X), len(X)))


def nearest(X, p, k=1):
    d = np.hypot(X[:, 0] - p[0], X[:, 1] - p[1])
    return np.argsort(d)[:k]


def path(pred, s, t):
    out = [t]
    while out[-1] != s and out[-1] >= 0:
        out.append(pred[out[-1]])
    return out[::-1]


def turning(pts, x1=(-110, -60), x2=(60, 110)):
    def ang(a, b):
        m = (pts[:, 0] >= a) & (pts[:, 0] <= b)
        if m.sum() < 3:
            return np.nan
        A = np.polyfit(pts[m, 0], pts[m, 1], 1)
        return np.arctan(A[0])
    return ang(*x1) - ang(*x2)          # positive = bent towards the mass (for b > 0)


# ---------------------------------------------------------------- continuum comparators
def speed_det(nval, two):
    return B37.rate(nval, nval, True) if two else nval


def ray_turning(b, two, L=120.0):
    """Shoot a ray with index n(x) = 1/speed(N(x)) from (-L, b) to (L, b); return its turning angle (same fit windows)."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    def nidx(x, y):
        N = 1 - 2 * GM / np.sqrt(x * x + y * y + RS * RS)
        return 1.0 / speed_det(N, two)
    def grad(x, y, h=1e-4):
        return ((nidx(x + h, y) - nidx(x - h, y)) / (2 * h), (nidx(x, y + h) - nidx(x, y - h)) / (2 * h))
    def shoot(a0, dense=False):
        def f(s, z):
            x, y, px, py = z
            n = nidx(x, y); gx, gy = grad(x, y)
            return [px / n, py / n, gx, gy]
        n0 = nidx(-L, b)
        ev = lambda s, z: z[0] - L
        ev.terminal = True
        sol = solve_ivp(f, (0, 4 * L), [-L, b, n0 * np.cos(a0), n0 * np.sin(a0)], events=ev, max_step=0.5,
                        rtol=1e-9, atol=1e-11, dense_output=dense)
        return sol
    def g(a0):
        ye = shoot(a0).y_events[0]
        return ye[0][1] - b if len(ye) else np.nan
    # primary image: the ray starts upward (away from the mass) and is bent back down to y = b
    grid = np.linspace(0.0, 0.35, 36)
    vals = [g(a0) for a0 in grid]
    k = next(i for i in range(len(grid) - 1) if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0)
    a = brentq(g, grid[k], grid[k + 1], xtol=1e-12)
    sol = shoot(a)
    pts = np.stack([sol.y[0], sol.y[1]], 1)
    return float(turning(pts)), float(a)


def straight_delay(b, two, L=120.0, nx=200001):
    x = np.linspace(-L, L, nx)
    N = 1 - 2 * GM / np.sqrt(x * x + b * b + RS * RS)
    f = 1 / speed_det(N, two) - 1 / speed_det(1.0, two)
    return float(np.trapezoid(f, x))


# ---------------------------------------------------------------- parts
def part_isotropy(stochastic):
    sec = 16; res = []
    for s in range(NREAL):
        X, pairs = graph(s)
        out = {}
        for two in (True, False):
            n = lapse(X, False)
            A = build(X, pairs, n, two, stochastic, 1000 + s)
            src = nearest(X, (0, 0))[0]
            T = dijkstra(A, indices=src)
            d = np.hypot(X[:, 0] - X[src, 0], X[:, 1] - X[src, 1])
            m = (d > 95) & (d < 105) & np.isfinite(T)
            ang = np.arctan2(X[m, 1] - X[src, 1], X[m, 0] - X[src, 0])
            k = ((ang + np.pi) / (2 * np.pi) * sec).astype(int) % sec
            v = d[m] / T[m]
            out["two" if two else "one"] = [float(v[k == q].mean()) for q in range(sec)]
        res.append(out)
    summ = {}
    for rule in ("two", "one"):
        V = np.array([r[rule] for r in res]).mean(0)
        summ[rule] = dict(sector_speed=V.tolist(), spread=float((V.max() - V.min()) / V.mean()), mean=float(V.mean()))
    return summ


def part_exponent(stochastic, nreal=6):
    Ns = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    out = {}
    for two in (True, False):
        V = []
        for Nv in Ns:
            vs = []
            for s in range(nreal):
                X, pairs = graph(s)
                n = np.full(len(X), Nv)
                A = build(X, pairs, n, two, stochastic, 2000 + s + int(100 * Nv))
                src = nearest(X, (0, 0))[0]
                T = dijkstra(A, indices=src)
                d = np.hypot(X[:, 0] - X[src, 0], X[:, 1] - X[src, 1])
                m = (d > 95) & (d < 105)
                vs.append(np.mean(d[m] / T[m]))
            V.append(np.mean(vs))
        p = np.polyfit(np.log(Ns), np.log(V), 1)[0]
        out["two" if two else "one"] = dict(N=list(Ns), v=[float(x) for x in V], p=float(p))
    return out


def part_lens(stochastic, bs=(20, 40, 60), L=120.0, nrecv=7):
    rows = {}
    per_real = []
    for s in range(NREAL):
        X, pairs = graph(s)
        r = np.hypot(X[:, 0], X[:, 1])
        keep = r >= R_CORE
        rec = {}
        for two in (True, False):
            key = "two" if two else "one"
            A1 = build(X, pairs, lapse(X, True), two, stochastic, 3000 + s, keep=keep)
            A0 = build(X, pairs, lapse(X, False), two, stochastic, 4000 + s, keep=keep)
            for b in bs:
                src = nearest(X, (-L, b))[0]
                tg = nearest(X, (L, b), nrecv)
                T1, pr1 = dijkstra(A1, indices=src, return_predecessors=True)
                T0 = dijkstra(A0, indices=src)
                delay = float(np.mean(T1[tg] - T0[tg]))
                flat_speed = float(np.mean(np.hypot(*(X[tg] - X[src]).T) / T0[tg]))
                turns, rmin = [], []
                for t in tg:
                    pp = X[path(pr1, src, t)]
                    turns.append(turning(pp)); rmin.append(np.hypot(pp[:, 0], pp[:, 1]).min())
                rec[f"{key}_{b}"] = dict(delay=delay, flat_speed=flat_speed, turn=float(np.nanmean(turns)),
                                         rmin=float(np.min(rmin)))
        per_real.append(rec)
    for key in per_real[0]:
        d = np.array([p[key]["delay"] for p in per_real]); tr = np.array([p[key]["turn"] for p in per_real])
        rows[key] = dict(delay=float(d.mean()), delay_se=float(d.std(ddof=1) / np.sqrt(len(d))),
                         turn=float(tr.mean()), turn_se=float(tr.std(ddof=1) / np.sqrt(len(tr))),
                         flat_speed=float(np.mean([p[key]["flat_speed"] for p in per_real])),
                         rmin=float(min(p[key]["rmin"] for p in per_real)),
                         halfA=dict(delay=float(d[:NREAL // 2].mean()), turn=float(tr[:NREAL // 2].mean())),
                         halfB=dict(delay=float(d[NREAL // 2:].mean()), turn=float(tr[NREAL // 2:].mean())))
    (OUT / f"lens_rows_{'stoch' if stochastic else 'det'}.json").write_text(json.dumps(rows, indent=1, default=float) + "\n")
    summ = {}
    for b in bs:
        t, o = rows[f"two_{b}"], rows[f"one_{b}"]
        ct, co = ray_turning(b, True)[0], ray_turning(b, False)[0]
        dt_c, do_c = straight_delay(b, True), straight_delay(b, False)
        # delays in flat travel-time units (as in 37B): delay x flat speed = equivalent extra length
        eq = lambda r, h=None: (r[h]["delay"] if h else r["delay"]) * r["flat_speed"]
        summ[b] = dict(delay_ratio=eq(t) / eq(o),
                       delay_ratio_halves=[eq(t, h) / eq(o, h) for h in ("halfA", "halfB")],
                       turn_ratio=t["turn"] / o["turn"],
                       turn_ratio_halves=[t[h]["turn"] / o[h]["turn"] for h in ("halfA", "halfB")],
                       continuum_turn_two=ct, continuum_turn_one=co, continuum_turn_ratio=ct / co,
                       turn_two_over_continuum=t["turn"] / ct, turn_one_over_continuum=o["turn"] / co,
                       straight_delay_ratio_continuum=(dt_c * speed_det(1.0, True)) / (do_c * speed_det(1.0, False)),
                       delay_two_equiv_length_over_straight_continuum=t["delay"] * t["flat_speed"] / (dt_c * speed_det(1.0, True)),
                       delay_one_equiv_length_over_straight_continuum=o["delay"] * o["flat_speed"] / (do_c * speed_det(1.0, False)))
    return dict(rows=rows, summary=summ)


def part_checks():
    out = {}
    # G5: Dijkstra on a chain reproduces the 37B rate sum along that chain
    x = np.arange(-120, 121, 1.0); b = 20.0
    X = np.stack([x, np.full_like(x, b)], 1)
    pairs = np.stack([np.arange(len(x) - 1), np.arange(1, len(x))], 1)
    for two in (True, False):
        n = lapse(X, True)
        A = build(X, pairs, n, two, False, 0)
        T = dijkstra(A, indices=0)[-1]
        out[f"G5_chain_{'two' if two else 'one'}"] = dict(dijkstra=float(T), rate_sum=B37.predicted_time(n, two),
                                                          rel=float(abs(T / B37.predicted_time(n, two) - 1)))
    # G6: mean sampled crossing time vs 1/rate
    g6 = {}
    for a in (0.5, 0.8, 1.0):
        smp = hop_samples(np.full(200000, a), np.full(200000, a), True, TAU, 7)
        g6[a] = dict(mean=float(smp.mean()), se=float(smp.std() / np.sqrt(len(smp))), inv_rate=float(1 / B37.rate(a, a, True)),
                     rel=float(smp.mean() * B37.rate(a, a, True) - 1))
    out["G6"] = g6
    return out


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.process_time()
    if part == "checks":
        res = part_checks()
    else:
        kind, mode = part.split("_")
        st = mode == "stoch"
        res = dict(isotropy=part_isotropy, exponent=part_exponent, lens=part_lens)[kind](st)
    res = dict(result=res, cpu_seconds=time.process_time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, indent=1, default=float)[:3000])
