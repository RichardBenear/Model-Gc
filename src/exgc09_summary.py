"""EX-GC09 summary: survival vs Omega, and the far G deficit against the lattice Green's function of a point sink.

Separate analysis path: reads only the JSON records; the Green's function is solved here (scipy sparse), not taken
from the run code.
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC09")
FIG = os.path.join(e.ROOT, "paper", "figures", "EX-GC09_token_lattice.png")


def green_profile(h, L, nb=40):
    """Deficit 1-G for a unit point sink at the centre cell, Dirichlet 0 just outside the cube (7-point stencil)."""
    n = int(round(2 * L / h)) + 1
    I = sp.identity(n, format="csr")
    d1 = sp.diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1], format="csr") / h ** 2
    A = sp.kron(sp.kron(d1, I), I) + sp.kron(sp.kron(I, d1), I) + sp.kron(sp.kron(I, I), d1)
    b = np.zeros(n ** 3)
    c = (n - 1) // 2
    b[(c * n + c) * n + c] = 1.0 / h ** 3          # -lap(phi) = delta  ->  phi = 1/(4 pi r) in free space
    phi = spl.spsolve((-A).tocsc(), b).reshape(n, n, n)
    ax = (np.arange(n) - c) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    edges = np.linspace(0, R.max(), nb + 1)
    idx = np.digitize(R.ravel(), edges) - 1
    out = []
    for k in range(nb):
        sel = idx == k
        if sel.sum() > 5:
            out.append((float(R.ravel()[sel].mean()), float(phi.ravel()[sel].mean())))
    return np.array(out)


def far_ratio(prof, M, gp, rmin, rmax):
    p = np.array(prof)
    sel = (p[:, 0] > rmin) & (p[:, 0] < rmax)
    pred = M * np.interp(p[sel, 0], gp[:, 0], gp[:, 1])
    return float(np.mean((1 - p[sel, 1]) / pred)), p[sel, 0], 1 - p[sel, 1], pred


def main():
    summ = dict(det={}, tok={})
    gps = {L: green_profile(0.5, L) for L in (6, 10)}
    for L in (6, 10):
        d = json.load(open(os.path.join(OUT, f"det_h0.5_L{L}.json")))
        M = d["rows"][-1]["M"]
        q, *_ = far_ratio(d["Gprofile"], M, gps[L], 3.5, L - 1)
        summ["det"][f"h0.5_L{L}"] = dict(M=M, far_ratio=q, T=d["rows"][-1]["t"])
    d = json.load(open(os.path.join(OUT, "det_h0.4_L10.json")))
    summ["det"]["h0.4_L10"] = dict(M=d["rows"][-1]["M"], T=d["rows"][-1]["t"])
    for Om in (100000, 1000, 300, 100, 30):
        d = json.load(open(os.path.join(OUT, f"tok_h0.5_O{Om}_s0.json")))
        rows = d["rows"]
        t = np.array([r["t"] for r in rows])
        M = np.array([r["M"] for r in rows])
        alive = rows[-1]["Wtot"] > 0
        late = t >= 10
        rec = dict(T=float(t[-1]), alive=bool(alive), clipped=rows[-1]["clipped"])
        if late.any():
            rec.update(M_mean=float(M[late].mean()), M_sd=float(M[late].std()), M_last=float(M[-1]))
            if alive:
                rec["far_ratio"] = far_ratio(d["Gprofile"], rec["M_mean"], gps[6], 3.5, 5)[0]
        else:
            rec.update(M_last=float(M[-1]))
        summ["tok"][str(Om)] = rec
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    det6 = json.load(open(os.path.join(OUT, "det_h0.5_L6.json")))
    ax[0].plot([r["t"] for r in det6["rows"]], [r["M"] for r in det6["rows"]], "k-", lw=2, label="deterministic")
    for Om, col in ((100000, "C4"), (1000, "C0"), (300, "C1"), (100, "C2"), (30, "C3")):
        d = json.load(open(os.path.join(OUT, f"tok_h0.5_O{Om}_s0.json")))
        ax[0].plot([r["t"] for r in d["rows"]], [r["M"] for r in d["rows"]], color=col, label=f"tokens, Ω = {Om:g}")
    ax[0].set_xlabel("time")
    ax[0].set_ylabel("G consumption M")
    ax[0].set_title("Gc1 particle on a 3D lattice (h = 0.5, box ±6)")
    ax[0].legend(fontsize=8)
    d10 = json.load(open(os.path.join(OUT, "det_h0.5_L10.json")))
    M10 = d10["rows"][-1]["M"]
    _, r, df, pr = far_ratio(d10["Gprofile"], M10, gps[10], 0.4, 9.5)
    ax[1].semilogy(r, df, "ko", ms=4, label="deterministic lattice, box ±10")
    ax[1].semilogy(r, pr, "k-", label="M × point-sink Green's function")
    for Om, col in ((1000, "C0"), (300, "C1")):
        d = json.load(open(os.path.join(OUT, f"tok_h0.5_O{Om}_s0.json")))
        Mm = summ["tok"][str(Om)]["M_mean"]
        _, r, df, pr = far_ratio(d["Gprofile"], Mm, gps[6], 0.4, 5.6)
        ax[1].semilogy(r, df, "o", color=col, ms=4, label=f"tokens Ω = {Om}, box ±6 (time average)")
        ax[1].semilogy(r, pr, "-", color=col, lw=0.8)
    ax[1].set_xlabel("distance r")
    ax[1].set_ylabel("G deficit 1 − G")
    ax[1].set_title("Far G field: sink of strength M")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130)


if __name__ == "__main__":
    main()
