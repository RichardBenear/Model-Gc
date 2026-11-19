"""EX-GC12 (exploratory): is a degree-1 frame texture (lattice 'Skyrmion') stable under the lift rule?

Frames on an n^3 lattice with the outer face fixed at the identity; each node also carries its SU(2) lift (one of the
48 elements of the binary octahedral group). Moves change one node's frame by a 90 or 120 degree step under
EX-GC08's move table ('lift' = three-node rule, 'flag' = pairwise control); the new lift is the one a short step away
from the old. Degree = signed count of Kuhn tetrahedra whose spherical image covers a generic point of S^3.
Our construction.
"""
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc08_belt as b  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC12")


def _lifts():
    Q = []
    for i in range(b.N):
        q = b.quat_short(i)
        if q is None:                                  # 180 degrees: axis from (m + I)/2
            M = (b.G[i].astype(float) + np.eye(3)) / 2
            ax = M[np.argmax(np.diag(M))]
            q = np.array([0.0, *(ax / np.linalg.norm(ax))])
        Q.append(q)
        Q.append(-q)
    return np.array(Q)


Q = _lifts()                                           # Q[2i], Q[2i+1] = the two lifts of frame i
MISM = np.array([[b.MUL[b.INV[a], c] != b.E_ID for c in range(b.N)] for a in range(b.N)])
NB = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
ID_LIFT = 2 * b.E_ID                                   # Q[2*E_ID] = +1


def hedgehog(n, R, centres=((0, 0, 0),), signs=(1,)):
    L = (n - 1) // 2
    ax = np.arange(n) - L
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    q = np.zeros((n, n, n, 4))
    q[..., 0] = 1
    for c, sg in zip(centres, signs):
        d = np.stack([X - c[0], Y - c[1], Z - c[2]], -1).astype(float)
        r = np.linalg.norm(d, axis=-1)
        F = np.where(r < R, np.pi * (1 - r / R), 0.0)
        rh = d / np.maximum(r, 1e-9)[..., None]
        qc = np.concatenate([np.cos(F)[..., None], sg * np.sin(F)[..., None] * rh], -1)
        # product (textures with disjoint supports just multiply)
        w1, v1 = q[..., 0], q[..., 1:]
        w2, v2 = qc[..., 0], qc[..., 1:]
        q = np.concatenate([(w1 * w2 - (v1 * v2).sum(-1))[..., None],
                            w1[..., None] * v2 + w2[..., None] * v1 + np.cross(v1, v2)], -1)
    return np.argmax(q @ Q.T, axis=-1)


def degree(lift, p):
    q = Q[lift]
    n = lift.shape[0]
    base = np.array(list(itertools.product(range(n - 1), repeat=3)))
    deg = 0
    for perm in itertools.permutations(range(3)):
        V = [np.zeros(3, int)]
        for a in perm:
            v = V[-1].copy()
            v[a] += 1
            V.append(v)
        Mx = np.stack([q[tuple((base + v).T)] for v in V], -1)
        det = np.linalg.det(Mx)
        ok = np.abs(det) > 1e-9
        c = np.linalg.solve(Mx[ok], np.broadcast_to(p, (int(ok.sum()), 4))[..., None])[..., 0]
        inside = (c > 0).all(-1)
        deg += int((np.sign(det[ok][inside]) * np.linalg.det(np.eye(3)[list(perm)])).sum())
    return deg


def min_dot(lift):
    q = Q[lift]
    out = 1.0
    for a in range(3):
        s = [slice(None)] * 3
        t = [slice(None)] * 3
        s[a], t[a] = slice(0, -1), slice(1, None)
        out = min(out, float((q[tuple(s)] * q[tuple(t)]).sum(-1).min()))
    return out


class Field:
    def __init__(self, lift, rule="lift", Tf=0.6, seed=0):
        self.S, self.A = b.move_table(("90", "120"), rule)
        self.S = np.array(self.S)
        self.lift = lift.copy()
        n = lift.shape[0]
        self.n = n
        ax = np.arange(n)
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        outer = (X == 0) | (X == n - 1) | (Y == 0) | (Y == n - 1) | (Z == 0) | (Z == n - 1)
        self.free = ~outer
        self.par = (X + Y + Z) % 2 == 0
        self.Tf = Tf
        self.rng = np.random.default_rng(seed)

    @property
    def f(self):
        return self.lift // 2

    def energy(self):
        f = self.f
        return int(MISM[f[:-1], f[1:]].sum() + MISM[f[:, :-1], f[:, 1:]].sum() + MISM[f[:, :, :-1], f[:, :, 1:]].sum())

    def sweep(self):
        n = self.n
        for colour in (True, False):
            f = self.f
            pad = np.pad(f, 1, mode="edge")
            nbr = np.stack([pad[1 + a:1 + a + n, 1 + c:1 + c + n, 1 + d:1 + d + n] for a, c, d in NB])
            idx = np.nonzero(self.free & (self.par == colour))
            g = f[idx]
            h = b.MUL[g, self.S[self.rng.integers(0, len(self.S), len(g))]]
            fn = nbr[:, idx[0], idx[1], idx[2]]
            ok = self.A[fn, g[None, :], h[None, :]].all(axis=0)
            dE = (MISM[h[None, :], fn].astype(int) - MISM[g[None, :], fn]).sum(axis=0)
            go = ok & ((dE <= 0) | (self.rng.random(len(g)) < np.exp(-np.maximum(dE, 0) / self.Tf)))
            old = self.lift[idx][go]
            hh = h[go]
            plus = 2 * hh
            dot = (Q[old] * Q[plus]).sum(-1)
            new = np.where(dot > 0, plus, plus + 1)
            sel = tuple(i[go] for i in idx)
            self.lift[sel] = new


def run(case, R=10, n=25, Tf=0.6, sweeps=10000, rule="lift", seed=0, every=500):
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(12345)
    pts = [p / np.linalg.norm(p) for p in rng.normal(size=(3, 4))]
    if case == "deg1":
        lift = hedgehog(n, R)
    elif case == "deg0":
        lift = np.full((n, n, n), ID_LIFT)
    elif case == "pair":                              # +1 and -1 textures side by side
        lift = hedgehog(n, 5.5, centres=((-6, 0, 0), (6, 0, 0)), signs=(1, -1))
    elif case == "two":                               # two +1 textures
        lift = hedgehog(n, 5.5, centres=((-6, 0, 0), (6, 0, 0)), signs=(1, 1))
    fd = Field(lift, rule=rule, Tf=Tf, seed=seed)
    hist = []
    t0 = time.time()
    for sw in range(sweeps + 1):
        if sw % every == 0:
            hist.append(dict(sweep=sw, E=fd.energy(), deg=[degree(fd.lift, p) for p in pts],
                             min_dot=round(min_dot(fd.lift), 3)))
        if sw < sweeps:
            fd.sweep()
    rec = dict(case=case, R=R, n=n, Tf=Tf, sweeps=sweeps, rule=rule, seed=seed, hist=hist, wall=time.time() - t0)
    json.dump(rec, open(os.path.join(OUT, f"{case}_{rule}_Tf{Tf:g}_s{seed}.json"), "w"), indent=1)
    h = hist
    print("RESULT", case, rule, Tf, seed, "E", h[0]["E"], "->", h[-1]["E"], "min", min(x["E"] for x in h),
          "deg", h[0]["deg"], "->", h[-1]["deg"], "degs seen", sorted({tuple(x["deg"]) for x in h}),
          "min_dot", min(x["min_dot"] for x in h), "wall", round(rec["wall"]), flush=True)
    return rec


if __name__ == "__main__":
    a = sys.argv
    run(a[1], Tf=float(a[2]), rule=a[3] if len(a) > 3 else "lift", seed=int(a[4]) if len(a) > 4 else 0,
        sweeps=int(a[5]) if len(a) > 5 else 10000)
