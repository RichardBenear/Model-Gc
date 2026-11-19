"""EX-GC13 (exploratory): Z2 loop class (pi_4(S^3)) of loops of SU(2) textures, by Pontryagin-Thom.

For a loop U_t (t in [0,1], U_1 = U_0) form V_t(x) = U_t(x) U_0(x)^-1. The preimage of a regular value p (near -1)
in the strip (0,1) x R^3 = R^4 is a set of closed curves. Each curve is framed by pulling back a basis of T_pS^3.
Score per component: class in pi_1(SO(4)) of (tangent, frame), read through the self-dual part, plus 1 (the round
circle with its bounding framing has class 1). The loop class is the sum mod 2.
Validation loops and the exchange test use compactly supported hedgehogs. Our construction of the test; the theorem
being tested is Finkelstein-Rubinstein (1968).
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC13")
TS = 8.0                                     # length scale of t in the R^4 metric


# ------------------------------------------------------------------ quaternions (w, x, y, z), last axis
def qmul(a, b):
    w1, x1, y1, z1 = np.moveaxis(a, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(b, -1, 0)
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], -1)


def qconj(a):
    return a * np.array([1, -1, -1, -1])


def rotz(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


# ------------------------------------------------------------------ textures
def hedgehog(x, R=2.0, sign=1):
    """Compactly supported degree-(sign) hedgehog: F = pi (1 - (r/R)^2)^3 for r < R."""
    r = np.linalg.norm(x, axis=-1)
    s = np.clip(r / R, 0, 1)
    F = np.pi * (1 - s ** 2) ** 3
    rh = x / np.maximum(r, 1e-12)[..., None]
    return np.concatenate([np.cos(F)[..., None], sign * np.sin(F)[..., None] * rh], -1)


def make_loop(kind, d=4.0, R=2.0):
    """Return U(t, x) for t in [0,1] (x: (..., 3))."""
    def rot_apply(th, y):
        return y @ rotz(th)                        # = R(th)^-1 y for row vectors

    if kind == "rot1":                              # 2 pi rotation of B = 1
        return lambda t, x: hedgehog(rot_apply(2 * np.pi * t, x), R)
    if kind == "rot2turns":                         # 4 pi rotation of B = 1
        return lambda t, x: hedgehog(rot_apply(4 * np.pi * t, x), R)
    if kind == "rotB2":                             # 2 pi rotation of B = 2 (pointwise square)
        def U(t, x):
            h = hedgehog(rot_apply(2 * np.pi * t, x), R)
            return qmul(h, h)
        return U
    if kind == "translate":                         # B = 1 carried around a circle, orientation fixed
        return lambda t, x: hedgehog(x - d * np.array([np.cos(2 * np.pi * t) - 1, np.sin(2 * np.pi * t), 0]), R)
    if kind in ("exchange", "exchange_rotA", "exchange_twice"):
        nt = 2.0 if kind == "exchange_twice" else 1.0

        def U(t, x):
            ph = np.pi * nt * t
            a = d * np.array([-np.cos(ph), -np.sin(ph), 0])
            bpos = -a
            ya = x - a
            if kind == "exchange_rotA":
                ya = rot_apply(2 * np.pi * t, ya)
            return qmul(hedgehog(ya, R), hedgehog(x - bpos, R))
        return U
    if kind == "rotpair":                           # rigid 2 pi rotation of the B = 2 pair about its midpoint
        def U(t, x):
            y = rot_apply(2 * np.pi * t, x)
            return qmul(hedgehog(y - np.array([-d, 0, 0]), R), hedgehog(y - np.array([d, 0, 0]), R))
        return U
    raise ValueError(kind)


# ------------------------------------------------------------------ preimage curves
class Problem:
    def __init__(self, U, p):
        self.U = U
        self.pc = qconj(p)

    def V(self, z):
        """z = (TS*t, x, y, z) -> V_t(x) (handles leading batch dims)."""
        t = z[..., 0] / TS
        x = z[..., 1:]
        tt = np.atleast_1d(t).ravel()
        xx = x.reshape(-1, 3)
        out = np.empty((len(xx), 4))
        for tv in np.unique(tt):
            sel = tt == tv
            out[sel] = qmul(self.U(tv, xx[sel]), qconj(self.U(0.0, xx[sel])))
        return out.reshape(z.shape[:-1] + (4,))

    def G(self, z):
        w = qmul(self.pc, self.V(z))
        return w[..., 1:], w[..., 0]

    def J(self, z, eps=1e-6):
        E = np.eye(4) * eps
        zz = np.concatenate([z[None] + E, z[None] - E])
        g, _ = self.G(zz)
        return ((g[:4] - g[4:]) / (2 * eps)).T            # 3 x 4

    def newton(self, z, fix_t=False, tol=1e-10, maxit=30):
        for _ in range(maxit):
            g, w = self.G(z)
            if np.linalg.norm(g) < tol:
                return z, w > 0
            Jm = self.J(z)
            if fix_t:
                dz = np.zeros(4)
                dz[1:] = np.linalg.lstsq(Jm[:, 1:], -g, rcond=None)[0]
            else:
                dz = np.linalg.pinv(Jm) @ (-g)
            z = z + dz
        g, w = self.G(z)
        return (z if np.linalg.norm(g) < 1e-7 else None), w > 0


def frame(Jm):
    """Tangent T (oriented so det[T, F] > 0) and frame F = J^+ (columns) at a curve point."""
    u, s, vt = np.linalg.svd(Jm)
    T = vt[-1]
    F = np.linalg.pinv(Jm)                               # 4 x 3
    if np.linalg.det(np.column_stack([T, F])) < 0:
        T = -T
    return T, F


def trace(pb, z0, ds=0.03, maxsteps=200000):
    T0, _ = frame(pb.J(z0))
    pts = [z0]
    z = z0
    T = T0
    for k in range(maxsteps):
        zp = z + ds * T
        # corrector on the hyperplane orthogonal to T through the predictor
        for _ in range(20):
            g, _ = pb.G(zp)
            Jm = pb.J(zp)
            A = np.vstack([Jm, T])
            rhs = np.concatenate([-g, [-(zp - z - ds * T) @ T]])
            dz = np.linalg.solve(A, rhs)
            zp = zp + dz
            if np.linalg.norm(dz) < 1e-11:
                break
        Tn, _ = frame(pb.J(zp))
        if Tn @ T < 0:
            raise RuntimeError("orientation flip on curve")
        z, T = zp, Tn
        pts.append(z)
        if k > 20 and np.linalg.norm(z - z0) < 1.5 * ds:
            return np.array(pts), True
    return np.array(pts), False


def selfdual(O):
    """SO(4) -> SO(3) via the action on self-dual 2-forms."""
    def W(a, b, c, e):
        m = np.zeros((4, 4))
        m[a, b], m[b, a] = 1, -1
        m[c, e], m[e, c] = 1, -1
        return m / 2.0
    om = [W(0, 1, 2, 3), W(0, 2, 3, 1), W(0, 3, 1, 2)]
    return np.array([[np.sum(om[i] * (O @ om[j] @ O.T)) * 2 for j in range(3)] for i in range(3)])


def rot_to_quat(Rm):
    w = np.sqrt(max(0.0, 1 + np.trace(Rm))) / 2
    x = np.sqrt(max(0.0, 1 + Rm[0, 0] - Rm[1, 1] - Rm[2, 2])) / 2
    y = np.sqrt(max(0.0, 1 - Rm[0, 0] + Rm[1, 1] - Rm[2, 2])) / 2
    zz = np.sqrt(max(0.0, 1 - Rm[0, 0] - Rm[1, 1] + Rm[2, 2])) / 2
    x = np.copysign(x, Rm[2, 1] - Rm[1, 2])
    y = np.copysign(y, Rm[0, 2] - Rm[2, 0])
    zz = np.copysign(zz, Rm[1, 0] - Rm[0, 1])
    q = np.array([w, x, y, zz])
    return q / np.linalg.norm(q)


def component_class(pb, pts):
    """pi_1(SO(4)) class of (T, F) along a closed curve (via the self-dual SO(3)), 0 or 1."""
    lift0 = None
    prev = None
    for z in pts:
        T, F = frame(pb.J(z))
        M = np.column_stack([T, F])
        Q_, Rr = np.linalg.qr(M)
        Q_ = Q_ * np.sign(np.diag(Rr))                   # Gram-Schmidt with positive diagonal
        if np.linalg.det(Q_) < 0:
            raise RuntimeError("frame orientation")
        q = rot_to_quat(selfdual(Q_))
        if prev is not None and q @ prev < 0:
            q = -q
        if lift0 is None:
            lift0 = q
        prev = q
    return int(prev @ lift0 < 0)


def find_curves(pb, box, tgrid=41, h=0.25, tol=0.5):
    curves = []
    ax = np.arange(-box, box + 1e-9, h)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    xs = np.stack([X.ravel(), Y.ravel(), Z.ravel()], -1)
    for t in np.linspace(0.04, 0.96, tgrid):
        g, w = pb.G(np.concatenate([np.full((len(xs), 1), TS * t), xs], -1))
        cand = xs[(np.linalg.norm(g, axis=-1) < tol) & (w > 0)]
        for x0 in cand:
            z, pos = pb.newton(np.concatenate([[TS * t], x0]), fix_t=False)
            if z is None or not pos or not (0 < z[0] < TS):
                continue
            if any(np.min(np.linalg.norm(c - z, axis=1)) < 0.2 for c in curves):
                continue
            pts, closed = trace(pb, z)
            if not closed:
                raise RuntimeError("curve did not close")
            curves.append(pts)
    return curves


def loop_class(kind, seed=0, d=4.0, R=2.0):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=4) * 0.5
    v[0] = -1.0
    p = v / np.linalg.norm(v)
    pb = Problem(make_loop(kind, d=d, R=R), p)
    box = d + R + 0.5 if kind not in ("rot1", "rot2turns", "rotB2") else R + 0.5
    t0 = time.time()
    curves = find_curves(pb, box)
    comps = [component_class(pb, c) for c in curves]
    total = int((sum(comps) + len(comps)) % 2)
    return dict(kind=kind, seed=seed, p=p.tolist(), n_components=len(curves), component_classes=comps,
                lengths=[len(c) for c in curves], loop_class=total, wall=time.time() - t0)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    kind = sys.argv[1]
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    r = loop_class(kind, seed)
    json.dump(r, open(os.path.join(OUT, f"loop_{kind}_s{seed}.json"), "w"), indent=1)
    print("RESULT", kind, seed, "components", r["n_components"], r["component_classes"], "class", r["loop_class"],
          "wall", round(r["wall"]), flush=True)
