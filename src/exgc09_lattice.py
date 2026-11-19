"""EX-GC09 (exploratory): Gc1 (G feeds X, k = 0.2) on a 3D cubic lattice graph, deterministic and as tokens.

Cells of spacing h; X, Y, W reflecting at the outer faces; G held at G_inf in a reservoir layer just outside the
faces (G conserved exactly in the bulk: only the core reaction G + W -> W + X removes it).
Token version: Omega tokens per unit concentration per unit volume, so a cell holds n = Omega h^3 c tokens.
Tau-leap with binomial hops/losses (counts never negative); production steps Poisson.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC09")
K = 0.2
D = np.array([1.0, e.DY, e.DW, 1.0])


def grid(h, L):
    n = int(round(2 * L / h)) + 1
    ax = (np.arange(n) - (n - 1) / 2) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    return ax, np.sqrt(X ** 2 + Y ** 2 + Z ** 2)


def init(h, L):
    z = np.load(os.path.join(e.ROOT, "results", "EX-GC01", "radial_ref.npz"))
    r0, u0 = z["r"], z["u"]
    ax, R = grid(h, L)
    return ax, R, np.array([np.interp(R, r0, u0[s]) for s in range(4)])


def lap(f, h, dirichlet=None):
    """7-point Laplacian; reflecting faces, or Dirichlet value 'dirichlet' just outside the faces."""
    if dirichlet is None:
        p = np.pad(f, 1, mode="edge")
    else:
        p = np.pad(f, 1, mode="constant", constant_values=dirichlet)
    return (p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1] + p[1:-1, :-2, 1:-1] + p[1:-1, 1:-1, 2:]
            + p[1:-1, 1:-1, :-2] - 6 * f) / h ** 2


def F(u, h):
    X, Y, W, G = u
    s = K * G * W
    return np.array([lap(X, h) + e.A - (1 + e.B) * X + X * X * Y + s,
                     e.DY * lap(Y, h) + e.B * X - X * X * Y,
                     e.DW * lap(W, h) - W + e.C * Y * W * W - W ** 3,
                     lap(G, h, e.GINF) - s])


def deterministic(h, L, T, dt):
    ax, R, u = init(h, L)
    rows = []
    for i in range(int(round(T / dt)) + 1):
        if i % int(round(1.0 / dt)) == 0:
            M = float(np.sum(K * u[3] * u[2]) * h ** 3)
            rows.append(dict(t=i * dt, M=M, W0=float(u[2].max()), finite=bool(np.isfinite(u).all())))
            if not rows[-1]["finite"]:
                break
        u = u + dt * F(u, h)
    return ax, R, u, rows


def gprofile(R, G, h, nb=40):
    edges = np.linspace(0, R.max(), nb + 1)
    idx = np.digitize(R.ravel(), edges) - 1
    out = []
    for b in range(nb):
        sel = idx == b
        if sel.sum() > 5:
            out.append((float(R.ravel()[sel].mean()), float(G.ravel()[sel].mean())))
    return out


# ---------------------------------------------------------------- tokens
def hop(n, rate_dt, rng, reflect=True, res=None):
    """Move tokens to the 6 neighbours. Each token leaves with prob 6*rate*dt (capped), direction uniform.
    reflect: hops through a face bounce back (stay). res: reservoir count per outer-face cell (Dirichlet G)."""
    p = min(6 * rate_dt, 0.9)
    leave = rng.binomial(n, p)
    out = n - leave
    rem = leave.copy()
    shares = []
    for k in range(5):
        x = rng.binomial(rem, 1.0 / (6 - k))
        shares.append(x)
        rem = rem - x
    shares.append(rem)
    new = out.copy()
    for k, (axis, sgn) in enumerate(((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1))):
        s = shares[k]
        moved = np.roll(s, sgn, axis=axis)
        # tokens that crossed the face: index 0 (sgn=+1 came from last) or last (sgn=-1 came from first)
        edge = [slice(None)] * 3
        src = [slice(None)] * 3
        if sgn == 1:
            edge[axis] = 0
            src[axis] = -1
        else:
            edge[axis] = -1
            src[axis] = 0
        moved[tuple(edge)] = 0
        if res is None:                      # reflecting: lost tokens stay in their cell
            new[tuple(src)] += s[tuple(src)]
        new += moved
    if res is not None:                      # Dirichlet reservoir: outflow lost, inflow Poisson from the reservoir
        for axis in range(3):
            for end in (0, -1):
                sl = [slice(None)] * 3
                sl[axis] = end
                new[tuple(sl)] += rng.poisson(res * rate_dt, size=new[tuple(sl)].shape)
    return new


def tokens(h, L, Omega, T, dt, seed=0, every=1.0, tavg=10.0):
    rng = np.random.default_rng(seed)
    ax, R, u = init(h, L)
    V = h ** 3
    OV = Omega * V
    n = np.rint(u * OV).astype(np.int64)
    n = np.maximum(n, 0)
    resG = e.GINF * OV
    rows = []
    nclip = 0
    Gacc = np.zeros(n[3].shape)
    nacc = 0
    t0 = time.time()
    for i in range(int(round(T / dt)) + 1):
        if i % int(round(every / dt)) == 0:
            Wc = n[2] / OV
            rows.append(dict(t=i * dt, M=float(np.sum(K * (n[3] / OV) * Wc) * V), Wtot=float(n[2].sum() / OV * V),
                             W0=float(Wc.max()), Gtot=int(n[3].sum()), clipped=nclip, wall=time.time() - t0))
            if i * dt >= tavg:
                Gacc += n[3] / OV
                nacc += 1
            if n[2].sum() == 0:
                break
        # Euler-consistent update: hops and reactions both use the state at the start of the step (hopping the
        # post-reaction counts instead gives an O(dt D/h^2) splitting bias, about 20% in M here). Removals are
        # drawn from the start-of-step counts and clipped against the post-hop counts.
        c = n / OV
        X0, Y0, W0n, G0 = n
        H = np.array([hop(n[0], D[0] / h ** 2 * dt, rng), hop(n[1], D[1] / h ** 2 * dt, rng),
                      hop(n[2], D[2] / h ** 2 * dt, rng), hop(n[3], D[3] / h ** 2 * dt, rng, res=resG)])
        prodX = rng.poisson(e.A * OV * dt, size=X0.shape)
        lossX = rng.binomial(X0, min(dt, 1.0))
        conv = rng.binomial(X0, min(e.B * dt, 1.0))                                  # B + X -> Y
        auto = np.minimum(rng.poisson(c[0] ** 2 * c[1] * OV * dt), Y0)               # 2X + Y -> 3X
        wl = rng.binomial(W0n, min(dt, 1.0))
        wg = rng.poisson(e.C * c[1] * c[2] ** 2 * OV * dt)                            # Y + 2W -> Y + 3W
        w3 = np.minimum(rng.poisson(c[2] ** 3 * OV * dt), W0n)                        # 3W -> 2W
        eat = np.minimum(rng.poisson(K * c[3] * c[2] * OV * dt), G0)                  # G + W -> W + X
        X = H[0] + prodX - lossX - conv + auto + eat
        Y = H[1] + conv - auto
        W = H[2] - wl - w3 + wg
        G = H[3] - eat
        clip = int(np.sum(np.minimum(X, 0)) + np.sum(np.minimum(Y, 0)) + np.sum(np.minimum(W, 0))
                   + np.sum(np.minimum(G, 0)))
        nclip += -clip
        n = np.maximum(np.array([X, Y, W, G]), 0)
    Gavg = Gacc / max(nacc, 1) if nacc else n[3] / OV
    return ax, R, n / OV, rows, Gavg


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    mode = sys.argv[1]
    if mode == "det":
        h, L, T = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
        dt = float(sys.argv[5]) if len(sys.argv) > 5 else 0.4 * h ** 2 / (6 * e.DY)
        ax, R, u, rows = deterministic(h, L, T, dt)
        prof = gprofile(R, u[3], h)
        json.dump(dict(h=h, L=L, dt=dt, rows=rows, Gprofile=prof), open(os.path.join(OUT, f"det_h{h:g}_L{L:g}.json"),
                                                                         "w"), indent=1)
        print("det h", h, "L", L, "final", rows[-1], flush=True)
    else:
        h, L, Om, T = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
        dt = float(sys.argv[7]) if len(sys.argv) > 7 else 0.002
        ax, R, c, rows, Gavg = tokens(h, L, Om, T, dt, seed)
        json.dump(dict(h=h, L=L, Omega=Om, dt=dt, seed=seed, rows=rows, Gprofile=gprofile(R, Gavg, h),
                       Gprofile_final=gprofile(R, c[3], h)),
                  open(os.path.join(OUT, f"tok_h{h:g}_O{Om:g}_s{seed}.json"), "w"), indent=1)
        print("tok h", h, "Omega", Om, "seed", seed, "final", rows[-1], flush=True)
