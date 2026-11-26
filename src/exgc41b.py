"""EX-GC41B: shared site clocks on a random pair network.

Each site has one unit-rate Poisson stream; its firing times are the stream's times divided by N_i (common random
numbers for flat and lens). A signal at x at time t crosses link (x,y) at the first coincidence after t (two-site:
one end fires while the other fired within tau before) or at x's first firing after t (one-site). Both are FIFO, so a
time-dependent Dijkstra gives the exact earliest arrival. Plan: spec/EX-GC41_plan.md. Our construction.
"""
import json, sys, time, os
from pathlib import Path
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc37b as B37
import exgc40b as B40

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC41"
TAU = B37.TAU
HALF, RC, RS, GM = 50.0, 1.5, 3.0, 1.0
XSRC, XREC = -45.0, 40.0
BS = (10, 20)
W = 8.0           # block length in unit time for lazily generated streams
NG, NH = 10, 21


@njit(cache=True)
def _block(seed, i, k):
    np.random.seed((seed * 1000003 + i * 7919 + k * 104729) % 2147483647)
    n = np.random.poisson(W)
    u = np.sort(np.random.random(n)) * W + k * W
    return u


@njit(cache=True)
def _last_before(seed, i, Ni, s):
    """Latest firing time of site i at or before s (real time); -inf if none in the 50 blocks searched back."""
    ku = int(np.floor(s * Ni / W))
    for k in range(ku, max(-1, ku - 50), -1):
        u = _block(seed, i, k)
        best = -np.inf
        for v in u:
            t = v / Ni
            if t <= s:
                best = t
        if best > -np.inf:
            return best
    return -np.inf


@njit(cache=True)
def _cross(seed, x, y, Nx, Ny, t, two, tau):
    """First crossing time of link x->y for a signal present at x from time t."""
    kx = int(np.floor(t * Nx / W)); ky = int(np.floor(t * Ny / W))
    if not two:
        k = kx
        while True:
            for v in _block(seed, x, k):
                tt = v / Nx
                if tt > t:
                    return tt
            k += 1
    # two-site: merge events of x and y after t; check the other end's latest firing within tau before
    lx = _last_before(seed, x, Nx, t); ly = _last_before(seed, y, Ny, t)
    ux = _block(seed, x, kx); uy = _block(seed, y, ky)
    ix = 0; iy = 0
    while True:
        while ix < len(ux) and ux[ix] / Nx <= t:
            ix += 1
        while iy < len(uy) and uy[iy] / Ny <= t:
            iy += 1
        if ix >= len(ux):
            kx += 1; ux = _block(seed, x, kx); ix = 0
            continue
        if iy >= len(uy):
            ky += 1; uy = _block(seed, y, ky); iy = 0
            continue
        tx = ux[ix] / Nx; ty = uy[iy] / Ny
        if tx <= ty:
            if tx - ly <= tau:
                return tx
            lx = tx; ix += 1
        else:
            if ty - lx <= tau:
                return ty
            ly = ty; iy += 1


@njit(cache=True)
def td_dijkstra(indptr, indices, Nsite, sources, seed, two, tau):
    n = len(Nsite)
    T = np.full(n, np.inf)
    done = np.zeros(n, np.bool_)
    # binary heap of (time, node)
    hk = np.empty(len(indices) + len(sources) + 1); hv = np.empty(len(indices) + len(sources) + 1, np.int64)
    hn = 0
    for s in sources:
        T[s] = 0.0
        hk[hn] = 0.0; hv[hn] = s; hn += 1
    while hn > 0:
        # pop min
        t = hk[0]; u = hv[0]
        hn -= 1
        hk[0] = hk[hn]; hv[0] = hv[hn]
        j = 0
        while True:
            l = 2 * j + 1; r = l + 1; m = j
            if l < hn and hk[l] < hk[m]:
                m = l
            if r < hn and hk[r] < hk[m]:
                m = r
            if m == j:
                break
            hk[j], hk[m] = hk[m], hk[j]; hv[j], hv[m] = hv[m], hv[j]; j = m
        if done[u]:
            continue
        done[u] = True
        for e in range(indptr[u], indptr[u + 1]):
            v = indices[e]
            if done[v]:
                continue
            tc = _cross(seed, u, v, Nsite[u], Nsite[v], t, two, tau)
            if tc < T[v]:
                T[v] = tc
                j = hn; hk[hn] = tc; hv[hn] = v; hn += 1
                while j > 0:
                    p = (j - 1) // 2
                    if hk[p] <= hk[j]:
                        break
                    hk[p], hk[j] = hk[j], hk[p]; hv[p], hv[j] = hv[j], hv[p]; j = p
    return T


def brute_force(indptr, indices, Nsite, sources, seed, two, tau):
    """H4: label-correcting relaxation over every link until nothing changes (no priority queue)."""
    n = len(Nsite); T = np.full(n, np.inf); T[sources] = 0.0
    changed = True
    while changed:
        changed = False
        for u in range(n):
            if not np.isfinite(T[u]):
                continue
            for e in range(indptr[u], indptr[u + 1]):
                v = indices[e]
                tc = _cross(seed, u, v, Nsite[u], Nsite[v], T[u], two, tau)
                if tc < T[v] - 1e-12:
                    T[v] = tc; changed = True
    return T


def graph(seed, half=HALF):
    rng = np.random.default_rng([20261126, 4102, seed])
    n = rng.poisson((2 * half) ** 2)
    X = (rng.random((n, 2)) * 2 - 1) * half
    pairs = cKDTree(X).query_pairs(RC, output_type="ndarray")
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    A = csr_matrix((np.ones(len(i)), (i, j)), shape=(n, n)); A.sort_indices()
    return X, pairs, A


def lapse(X, gm):
    return 1 - 2 * gm / np.sqrt(X[:, 0] ** 2 + X[:, 1] ** 2 + RS * RS)


def det_T(X, pairs, n, two, src):
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    w = 1.0 / B37.rate(n[i], n[j], True) if two else 1.0 / n[i]
    return dijkstra(csr_matrix((w, (i, j)), shape=(len(X), len(X))), indices=src, min_only=True)


def ind_T(X, pairs, n, two, src, seed):
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    w = B40.hop_samples(n[i].astype(float), n[j].astype(float), two, TAU, seed)
    return dijkstra(csr_matrix((w, (i, j)), shape=(len(X), len(X))), indices=src, min_only=True)


def readout(X, T0, T1, rec, band):
    bb = band[np.isfinite(T0[band])]
    vflat = 1.0 / np.polyfit(X[bb, 0], T0[bb], 1)[0]
    ok = rec[np.isfinite(T0[rec]) & np.isfinite(T1[rec])]
    y = X[ok, 1]; D = (T1[ok] - T0[ok]) * vflat
    return dict(vflat=float(vflat), **{str(b): float(D[np.abs(y - b) <= 5].mean()) for b in BS})


def one_graph(g):
    X, pairs, A = graph(g)
    src = np.where(X[:, 0] < XSRC)[0]
    rec = np.where((X[:, 0] >= XREC - 2) & (X[:, 0] <= XREC + 2))[0]
    band = np.where((np.abs(X[:, 1]) < 30) & (X[:, 0] > -35) & (X[:, 0] < 45))[0]
    n0 = np.ones(len(X)); n1 = lapse(X, GM)
    out = dict(det={}, ind=[], shared=[], median={})
    for two in (True, False):
        r = "two" if two else "one"
        out["det"][r] = readout(X, det_T(X, pairs, n0, two, src), det_T(X, pairs, n1, two, src), rec, band)
        S0, S1 = [], []
        for h in range(NH):
            seed = 1000 * g + h + 1
            if h < 5:
                out["ind"].append(dict(rule=r, h=h, **readout(X, ind_T(X, pairs, n0, two, src, 7 * seed),
                                                            ind_T(X, pairs, n1, two, src, 7 * seed + 3), rec, band)))
            t0 = td_dijkstra(A.indptr, A.indices, n0, src, seed, two, TAU)
            t1 = td_dijkstra(A.indptr, A.indices, n1, src, seed, two, TAU)
            S0.append(t0); S1.append(t1)
            out["shared"].append(dict(rule=r, h=h, **readout(X, t0, t1, rec, band)))
        out["median"][r] = readout(X, np.median(np.array(S0), 0), np.median(np.array(S1), 0), rec, band)
    return out


def checks():
    res = {}
    # H3: mean single-hop crossing time from a random start, shared-clock code vs the 40B sampler
    for a in (0.5, 1.0):
        ts = []
        for s in range(20000):
            t0 = 37.0 + (s % 97) * 1.37
            ts.append(_cross(s + 1, 0, 1, a, a, t0, True, TAU) - t0)
        smp = B40.hop_samples(np.full(200000, a), np.full(200000, a), True, TAU, 11)
        res[f"H3_{a}"] = dict(shared=float(np.mean(ts)), shared_se=float(np.std(ts) / np.sqrt(len(ts))),
                              sampler=float(smp.mean()), rel=float(np.mean(ts) / smp.mean() - 1))
    # H4: time-dependent Dijkstra equals brute-force relaxation on a 200-site graph
    rng = np.random.default_rng(5); X = rng.random((200, 2)) * 14 - 7
    pairs = cKDTree(X).query_pairs(RC, output_type="ndarray")
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    A = csr_matrix((np.ones(len(i)), (i, j)), shape=(200, 200)); A.sort_indices()
    n = lapse(X, 0.5); src = np.where(X[:, 0] < -6)[0]
    dev = 0.0
    for h in range(5):
        for two in (True, False):
            a = td_dijkstra(A.indptr, A.indices, n, src, 50 + h, two, TAU)
            b = brute_force(A.indptr, A.indices, n, src, 50 + h, two, TAU)
            fin = np.isfinite(a)
            assert np.array_equal(fin, np.isfinite(b))
            dev = max(dev, float(np.abs(a[fin] - b[fin]).max()))
    res["H4_max_dev"] = dev
    return res


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if part == "checks":
        res = checks()
    elif part == "run":
        from multiprocessing import Pool
        with Pool(int(os.environ.get("MODELG_WORKERS", "4"))) as p:
            per = p.map(one_graph, list(range(NG)))
        res = dict(per_graph=per)
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"41b_{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, default=float)[:1500])
