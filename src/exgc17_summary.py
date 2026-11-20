"""EX-GC17 summary: compare measured drag (gamma) and acceleration response (kappa) with a linear mean-field path.

Separate path: continuous-time mean token field on the same lattice, generator A(u) = A0 + u A1 (hop rates 1 along
x, y and 1 +- u/2 along +-z), Dirichlet reservoir inflow b(u) = b0 + u b1, fixed sinks m0 at the consumer sites.
To first order in u the field responds as n1(w) = (i w - A0)^-1 (A1 n0 + b1). The force is m0 times the net direction of
arrivals, linearised in u and n. gamma = Re F1, kappa = -Im F1 / w. Shares only make_units with the simulator.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc16_recoil as r16  # noqa: E402  (geometry only)

OUT = os.path.join(r16.ROOT, "results", "EX-GC17")
DIRS = ((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1))
RATE0 = np.ones(6)
RATE1 = np.array([0, 0, 0, 0, 0.5, -0.5])


def operators(L, Om):
    N = 2 * L + 1
    idx = np.arange(N ** 3).reshape(N, N, N)
    mats = [[], []]
    b = [np.zeros(N ** 3), np.zeros(N ** 3)]
    for k, (a, s) in enumerate(DIRS):
        src = [slice(None)] * 3
        dst = [slice(None)] * 3
        if s == 1:
            src[a], dst[a] = slice(0, -1), slice(1, None)
        else:
            src[a], dst[a] = slice(1, None), slice(0, -1)
        rows, cols = idx[tuple(dst)].ravel(), idx[tuple(src)].ravel()
        face = [slice(None)] * 3
        face[a] = 0 if s == 1 else N - 1
        fi = idx[tuple(face)].ravel()
        for j, rate in enumerate((RATE0[k], RATE1[k])):
            if rate:
                mats[j].append(sp.csr_matrix((np.full(rows.size, rate), (rows, cols)), shape=(N ** 3, N ** 3)))
                b[j][fi] += rate * Om
    A0 = sum(mats[0]) - 6 * sp.identity(N ** 3, format="csr")
    A1 = sum(mats[1])
    return idx, A0.tocsc(), A1.tocsc(), b[0], b[1]


def force(units, n, u):
    """z force for field n (flat-indexed via 3D array) and drift u."""
    F = 0.0
    for un in units:
        num, tot = 0.0, 0.0
        for s in un["sites"]:
            for k, (a, sg) in enumerate(DIRS):
                src = np.array(s)
                src[a] -= sg
                arr = (RATE0[k] + u * RATE1[k]) * n[tuple(src)]
                tot += arr
                if a == 2:
                    num += sg * arr
        F += un["rate"] * num / tot
    return F


def predict(kind, m0, d, L, Om, ws):
    units, B = r16.make_units(kind, L, m0, d)
    N = 2 * L + 1
    idx, A0, A1, b0, b1 = operators(L, Om)
    sink = np.zeros(N ** 3)
    for un in units:
        for s in un["sites"]:
            sink[idx[tuple(s)]] += un["rate"] / len(un["sites"])
    n0 = spla.spsolve(A0, -(b0 - sink))                        # 0 = A0 n0 + b0 - sink
    rhs = A1 @ n0 + b1
    eps = 1e-6
    n0g = n0.reshape(N, N, N)
    Fu = (force(units, n0g, eps) - force(units, n0g, -eps)) / (2 * eps)
    out = {}
    for w in ws:
        n1 = spla.spsolve((1j * w * sp.identity(N ** 3) - A0).tocsc(), rhs).reshape(N, N, N)
        dF = (force(units, n0g + eps * n1.real, 0) - force(units, n0g - eps * n1.real, 0)) / (2 * eps) \
            + 1j * (force(units, n0g + eps * n1.imag, 0) - force(units, n0g - eps * n1.imag, 0)) / (2 * eps)
        F1 = Fu + dF
        out[w] = dict(gamma_per_B=float(F1.real / B), kappa_per_B=float(-F1.imag / w / B) if w > 0 else None)
    return out


def main():
    files = sorted(glob.glob(os.path.join(OUT, "*_s*.json")))
    groups = defaultdict(list)
    for f in files:
        rec = json.load(open(f))
        groups[(rec["kind"], rec["m0"], rec["u"], rec["w"], rec["d"], rec["L"])].append(rec)
    rows = []
    cache = {}
    for (kind, m0, u, w, d, L), recs in sorted(groups.items()):
        key = (kind, m0, d, L)
        if key not in cache:
            cache[key] = predict(kind, m0, d, L, recs[0]["Om"], [0.0, 0.2, 1.0])
        p = cache[key].get(w) or predict(kind, m0, d, L, recs[0]["Om"], [w])[w]
        g = np.array([x["gamma_per_B"] for x in recs])
        row = dict(kind=kind, m0=m0, u=u, w=w, seeds=len(recs), gamma=float(g.mean()),
                   gamma_sem=float(g.std(ddof=1) / np.sqrt(len(g))) if len(g) > 1 else None,
                   gamma_pred=p["gamma_per_B"], gamma_weak=m0 / 6)
        line = (f"{kind:12s} m0={m0:<4g} u={u:<5g} w={w:<4g} n={len(recs)}  gamma/B={row['gamma']:8.3f}"
                f" (pred {p['gamma_per_B']:8.3f}, weak m0/6 {m0 / 6:7.3f}, ratio {row['gamma'] / (m0 / 6):.4f})")
        if w > 0:
            k = np.array([x["kappa_per_B"] for x in recs])
            row.update(kappa=float(k.mean()), kappa_sem=float(k.std(ddof=1) / np.sqrt(len(k))) if len(k) > 1 else None,
                       kappa_pred=p["kappa_per_B"])
            line += f"  kappa/B={row['kappa']:7.3f}" + (f" ±{row['kappa_sem']:.3f}" if row["kappa_sem"] else "") \
                + f" (pred {p['kappa_per_B']:7.3f})"
        print(line)
        rows.append(row)
    json.dump(rows, open(os.path.join(OUT, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
