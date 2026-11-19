"""EX-GC10 (exploratory): route 3, stage 3 - a Gc particle tethered by its own G supply.

Frames (one of the 24 cube rotations) sit on the EX-GC09 lattice nodes. Frame moves follow EX-GC08's three-node lift
rule (steps 90 + 120 degrees) and are attempted at a rate set by the G traffic into the node. A G hop across an edge
whose two frames differ succeeds with probability 1 - eps (twisted ribbons impede supply).
Chemistry: deterministic (explicit Euler, as EX-GC09 'det') or tokens (as EX-GC09 'tok', Euler-consistent update).
Our construction, not LaViolette's.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc08_belt as b  # noqa: E402
import exgc09_lattice as g9  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC10")
K = g9.K
D = g9.D
MISM = np.array([[b.MUL[b.INV[a], c] != b.E_ID for c in range(b.N)] for a in range(b.N)])
NB = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


class Frames:
    def __init__(self, n, turns, rule="lift", core=2, Tf=0.5, seed=0):
        self.S, self.A = b.move_table(("90", "120"), rule)
        self.S = np.array(self.S)
        self.n = n
        L = (n - 1) // 2
        self.L = L
        ax = np.arange(n) - L
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        layer = np.maximum(np.maximum(abs(X), abs(Y)), abs(Z))
        q = np.zeros_like(layer)
        for k in range(core + 1, L):
            q[layer == k] = int(round(4 * turns * (L - k) / (L - core - 1)))
        self.f = np.vectorize(b.rot_z)(q)
        self.f[layer <= core] = b.E_ID
        self.f[layer == L] = b.E_ID
        self.free = (layer > core) & (layer < L)
        self.par = ((X + Y + Z) % 2 == 0)
        self.Tf = Tf
        self.rng = np.random.default_rng(seed)

    def edge_mism(self):
        """Mismatch on the +x, +y, +z edges: arrays of shape (n-1,n,n), (n,n-1,n), (n,n,n-1)."""
        f = self.f
        return (MISM[f[:-1], f[1:]], MISM[f[:, :-1], f[:, 1:]], MISM[f[:, :, :-1], f[:, :, 1:]])

    def energy(self):
        return int(sum(m.sum() for m in self.edge_mism()))

    def by_layer(self):
        """Mismatched edges counted at the inner node's cube layer (|x|inf in lattice units)."""
        n = self.n
        ax = np.arange(n) - self.L
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        layer = np.maximum(np.maximum(abs(X), abs(Y)), abs(Z))
        out = np.zeros(self.L + 1, int)
        mx, my, mz = self.edge_mism()
        for m, lay in ((mx, np.minimum(layer[:-1], layer[1:])), (my, np.minimum(layer[:, :-1], layer[:, 1:])),
                       (mz, np.minimum(layer[:, :, :-1], layer[:, :, 1:]))):
            np.add.at(out, lay[m], 1)
        return out.tolist()

    def z2(self):
        return b.z2_lift(self.f, self.L)

    def update(self, p_attempt):
        """One checkerboard sweep; node attempts with probability p_attempt (array)."""
        f = self.f
        n = self.n
        acc = 0
        for colour in (True, False):
            pad = np.pad(f, 1, mode="edge")           # outer layer is fixed and never free, so padding is unused
            nbr = np.stack([pad[1 + a:1 + a + n, 1 + bb:1 + bb + n, 1 + c:1 + c + n] for a, bb, c in NB])
            sel = self.free & (self.par == colour) & (self.rng.random(f.shape) < p_attempt)
            idx = np.nonzero(sel)
            g = f[idx]
            h = b.MUL[g, self.S[self.rng.integers(0, len(self.S), len(g))]]
            fn = nbr[:, idx[0], idx[1], idx[2]]                  # (6, m)
            ok = self.A[fn, g[None, :], h[None, :]].all(axis=0)
            dE = (MISM[h[None, :], fn].astype(int) - MISM[g[None, :], fn]).sum(axis=0)
            u = self.rng.random(len(g))
            go = ok & ((dE <= 0) | (u < np.exp(-np.maximum(dE, 0) / self.Tf)))
            f[idx[0][go], idx[1][go], idx[2][go]] = h[go]
            acc += int(go.sum())
        return acc


def wlap(Gf, h, w, ginf):
    """Weighted Laplacian for G: edge conductance w (1 - eps*mismatch) inside; reservoir edges weight 1."""
    wx, wy, wz = w
    p = np.pad(Gf, 1, mode="constant", constant_values=ginf)
    out = np.zeros_like(Gf)
    # interior edges
    dx = Gf[1:] - Gf[:-1]
    out[:-1] += wx * dx
    out[1:] -= wx * dx
    dy = Gf[:, 1:] - Gf[:, :-1]
    out[:, :-1] += wy * dy
    out[:, 1:] -= wy * dy
    dz = Gf[:, :, 1:] - Gf[:, :, :-1]
    out[:, :, :-1] += wz * dz
    out[:, :, 1:] -= wz * dz
    # reservoir edges (the outer frames are fixed at identity, as is the reservoir: weight 1)
    out[0] += p[0, 1:-1, 1:-1] - Gf[0]
    out[-1] += p[-1, 1:-1, 1:-1] - Gf[-1]
    out[:, 0] += p[1:-1, 0, 1:-1] - Gf[:, 0]
    out[:, -1] += p[1:-1, -1, 1:-1] - Gf[:, -1]
    out[:, :, 0] += p[1:-1, 1:-1, 0] - Gf[:, :, 0]
    out[:, :, -1] += p[1:-1, 1:-1, -1] - Gf[:, :, -1]
    return out / h ** 2


def traffic(Gf, w):
    """G traffic into each node relative to the untwisted far field (G = 1, all 6 edges weight 1)."""
    n = Gf.shape[0]
    p = np.pad(Gf, 1, mode="constant", constant_values=e.GINF)
    wx, wy, wz = w
    t = np.zeros_like(Gf)
    t[1:] += wx * Gf[:-1]
    t[:-1] += wx * Gf[1:]
    t[:, 1:] += wy * Gf[:, :-1]
    t[:, :-1] += wy * Gf[:, 1:]
    t[:, :, 1:] += wz * Gf[:, :, :-1]
    t[:, :, :-1] += wz * Gf[:, :, 1:]
    t[0] += p[0, 1:-1, 1:-1]
    t[-1] += p[-1, 1:-1, 1:-1]
    t[:, 0] += p[1:-1, 0, 1:-1]
    t[:, -1] += p[1:-1, -1, 1:-1]
    t[:, :, 0] += p[1:-1, 1:-1, 0]
    t[:, :, -1] += p[1:-1, 1:-1, -1]
    assert t.shape == (n, n, n)
    return t / 6.0


def settled(h, L):
    """Deterministic settled particle on the lattice (cached)."""
    path = os.path.join(OUT, f"settled_h{h:g}_L{L:g}.npy")
    if os.path.exists(path):
        return np.load(path)
    ax, R, u, rows = g9.deterministic(h, L, 30.0, 0.0015)
    np.save(path, u)
    return u


def run(mode, turns, eps, rule="lift", Tf=0.5, nu=1.0, T=40.0, Omega=3000.0, h=0.5, L=6.0, seed=0, kf=10,
        dt=0.0015, every=1.0):
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(10_000 + seed)
    u = settled(h, L).copy()
    n = u.shape[1]
    fr = Frames(n, turns, rule=rule, Tf=Tf, seed=seed)
    V = h ** 3
    OV = Omega * V
    if mode == "tok":
        nt = np.maximum(np.rint(u * OV).astype(np.int64), 0)
    rows = []
    t0 = time.time()
    nsteps = int(round(T / dt))
    w = tuple(1.0 - eps * m for m in fr.edge_mism())
    acc = 0
    for i in range(nsteps + 1):
        c = u if mode == "det" else nt / OV
        if i % int(round(every / dt)) == 0:
            wsum = c[2].sum()
            cen = [float((c[2] * ((np.arange(n) - (n - 1) / 2) * h).reshape(sh)).sum() / max(wsum, 1e-12))
                   for sh in ((-1, 1, 1), (1, -1, 1), (1, 1, -1))]
            rows.append(dict(t=round(i * dt, 6), M=float(np.sum(K * c[3] * c[2]) * V), Wtot=float(c[2].sum() * V),
                             Wcen=[round(x, 4) for x in cen], Wmax=float(c[2].max()),
                             E=fr.energy(), z2=fr.z2(), acc=acc, wall=round(time.time() - t0, 1)))
            if c[2].sum() * V < 1e-3:
                break
        if i == nsteps:
            break
        if i % kf == 0:                                  # frame sweep, attempts driven by local G traffic
            acc += fr.update(np.minimum(1.0, nu * traffic(c[3], w)))
            w = tuple(1.0 - eps * m for m in fr.edge_mism())
        if mode == "det":
            X, Y, W, Gf = u
            s = K * Gf * W
            u = u + dt * np.array([g9.lap(X, h) + e.A - (1 + e.B) * X + X * X * Y + s,
                                   e.DY * g9.lap(Y, h) + e.B * X - X * X * Y,
                                   e.DW * g9.lap(W, h) - W + e.C * Y * W * W - W ** 3,
                                   wlap(Gf, h, w, e.GINF) - s])
        else:
            nt = tok_step(nt, c, OV, h, dt, w, eps, rng)
    tag = f"{mode}_t{turns}_e{eps:g}_{rule}_Tf{Tf:g}_s{seed}" + (f"_O{Omega:g}" if mode == "tok" else "")
    rec = dict(mode=mode, turns=turns, eps=eps, rule=rule, Tf=Tf, nu=nu, T=T, Omega=Omega if mode == "tok" else None,
               h=h, L=L, seed=seed, kf=kf, dt=dt, rows=rows, final_frames_E=fr.energy(),
               E_by_layer=fr.by_layer())
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1, default=lambda o: o.item())
    print("RESULT", tag, "E", rows[0]["E"], "->", rows[-1]["E"], "z2", rows[0]["z2"], "->", rows[-1]["z2"],
          "M", round(rows[0]["M"], 3), "->", round(rows[-1]["M"], 3), "wall", round(time.time() - t0), flush=True)
    return rec


def hop_w(nG, rate_dt, rng, w, res):
    """G hops with edge pass probability w (blocked tokens stay); reservoir inflow through unit-weight faces."""
    p = min(6 * rate_dt, 0.9)
    leave = rng.binomial(nG, p)
    new = nG - leave
    rem = leave.copy()
    shares = []
    for k in range(5):
        x = rng.binomial(rem, 1.0 / (6 - k))
        shares.append(x)
        rem = rem - x
    shares.append(rem)
    wx, wy, wz = w
    for k, (axis, sgn) in enumerate(((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1))):
        s = shares[k]
        we = (wx, wy, wz)[axis]
        # source slices: cells that have a neighbour in direction sgn along axis
        src = [slice(None)] * 3
        dst = [slice(None)] * 3
        face = [slice(None)] * 3
        if sgn == 1:
            src[axis], dst[axis], face[axis] = slice(0, -1), slice(1, None), -1
        else:
            src[axis], dst[axis], face[axis] = slice(1, None), slice(0, -1), 0
        s_in = s[tuple(src)]
        passed = rng.binomial(s_in, we) if np.any(we < 1) else s_in
        new[tuple(dst)] += passed
        new[tuple(src)] += s_in - passed                 # blocked: stay
        # tokens leaving through the outer face go to the reservoir (lost)
    for axis in range(3):
        for end in (0, -1):
            sl = [slice(None)] * 3
            sl[axis] = end
            new[tuple(sl)] += rng.poisson(res * rate_dt, size=new[tuple(sl)].shape)
    return new


def tok_step(n, c, OV, h, dt, w, eps, rng):
    X0, Y0, W0n, G0 = n
    H = np.array([g9.hop(n[0], D[0] / h ** 2 * dt, rng), g9.hop(n[1], D[1] / h ** 2 * dt, rng),
                  g9.hop(n[2], D[2] / h ** 2 * dt, rng), hop_w(n[3], D[3] / h ** 2 * dt, rng, w, e.GINF * OV)])
    prodX = rng.poisson(e.A * OV * dt, size=X0.shape)
    lossX = rng.binomial(X0, min(dt, 1.0))
    conv = rng.binomial(X0, min(e.B * dt, 1.0))
    auto = np.minimum(rng.poisson(c[0] ** 2 * c[1] * OV * dt), Y0)
    wl = rng.binomial(W0n, min(dt, 1.0))
    wg = rng.poisson(e.C * c[1] * c[2] ** 2 * OV * dt)
    w3 = np.minimum(rng.poisson(c[2] ** 3 * OV * dt), W0n)
    eat = np.minimum(rng.poisson(K * c[3] * c[2] * OV * dt), G0)
    return np.maximum(np.array([H[0] + prodX - lossX - conv + auto + eat, H[1] + conv - auto,
                                H[2] - wl - w3 + wg, H[3] - eat]), 0)


def frames_only(turns, Tf, sweeps, rule="lift", n=25, seed=0):
    fr = Frames(n, turns, rule=rule, Tf=Tf, seed=seed)
    hist = [(0, fr.energy(), fr.z2())]
    t0 = time.time()
    for sw in range(1, sweeps + 1):
        fr.update(np.ones((n, n, n)))
        if sw % 100 == 0:
            hist.append((sw, fr.energy(), fr.z2()))
    rec = dict(turns=turns, Tf=Tf, sweeps=sweeps, rule=rule, n=n, seed=seed, hist=hist, wall=time.time() - t0)
    json.dump(rec, open(os.path.join(OUT, f"frames_t{turns}_Tf{Tf:g}_{rule}_s{seed}.json"), "w"), indent=1,
              default=lambda o: o.item())
    print("FRAMES", turns, Tf, rule, seed, hist[0], "->", hist[-1], "min E", min(x[1] for x in hist),
          "wall", round(rec["wall"]), flush=True)
    return rec


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "frames":
        frames_only(int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]),
                    rule=sys.argv[5] if len(sys.argv) > 5 else "lift", seed=int(sys.argv[6]) if len(sys.argv) > 6
                    else 0)
    else:
        mode, turns, eps = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
        kw = {}
        for a in sys.argv[4:]:
            k, v = a.split("=")
            kw[k] = v if k == "rule" else float(v)
        if "seed" in kw:
            kw["seed"] = int(kw["seed"])
        if "kf" in kw:
            kw["kf"] = int(kw["kf"])
        run(mode, turns, eps, **kw)
