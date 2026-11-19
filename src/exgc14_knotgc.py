"""EX-GC14 (exploratory): the lattice knot (EX-GC12) as a Gc particle.

Deterministic Gc1 chemistry (EX-GC09 lattice, h = 0.5) plus frames with lifts (EX-GC12, lift rule).
Coupling (our design):
  (A) W autocatalysis c Y W^2 is multiplied by g = clip(rho_s / rho_ref, 0, 1), rho = (1 - Re q)/2 the knot density,
      rho_s its Gaussian smoothing (width sig, length units);
  (B) frame moves feel an extra energy -lam * sum W rho (knot drawn to the W core).
Optional external gradient: reservoir G = G_inf + E z.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc08_belt as b  # noqa: E402
import exgc09_lattice as g9  # noqa: E402
import exgc10_tether as t10  # noqa: E402
import exgc12_soliton as s12  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC14")
K = g9.K
RHO = (1 - s12.Q[:, 0]) / 2                     # knot density per lift


def lap_bc(f, h, bc):
    p = np.pad(f, 1, mode="constant")
    p[0, 1:-1, 1:-1], p[-1, 1:-1, 1:-1] = bc[0], bc[1]
    p[1:-1, 0, 1:-1], p[1:-1, -1, 1:-1] = bc[2], bc[3]
    p[1:-1, 1:-1, 0], p[1:-1, 1:-1, -1] = bc[4], bc[5]
    return (p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1] + p[1:-1, :-2, 1:-1] + p[1:-1, 1:-1, 2:]
            + p[1:-1, 1:-1, :-2] - 6 * f) / h ** 2


class KnotField(s12.Field):
    def __init__(self, lift, lam, **kw):
        super().__init__(lift, **kw)
        self.lam = lam
        self.Wf = np.zeros(lift.shape)

    def sweep(self):
        n = self.n
        for colour in (True, False):
            f = self.f
            pad = np.pad(f, 1, mode="edge")
            nbr = np.stack([pad[1 + a:1 + a + n, 1 + c:1 + c + n, 1 + d:1 + d + n] for a, c, d in s12.NB])
            idx = np.nonzero(self.free & (self.par == colour))
            g = f[idx]
            h = b.MUL[g, self.S[self.rng.integers(0, len(self.S), len(g))]]
            fn = nbr[:, idx[0], idx[1], idx[2]]
            ok = self.A[fn, g[None, :], h[None, :]].all(axis=0)
            old = self.lift[idx]
            plus = 2 * h
            new = np.where((s12.Q[old] * s12.Q[plus]).sum(-1) > 0, plus, plus + 1)
            dE = (s12.MISM[h[None, :], fn].astype(int) - s12.MISM[g[None, :], fn]).sum(axis=0).astype(float)
            dE -= self.lam * self.Wf[idx] * (RHO[new] - RHO[old])
            go = ok & ((dE <= 0) | (self.rng.random(len(g)) < np.exp(-np.maximum(dE, 0) / self.Tf)))
            sel = tuple(i[go] for i in idx)
            self.lift[sel] = new[go]


def centre(w, ax):
    s = w.sum()
    if s <= 0:
        return [float("nan")] * 3
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    return [float((w * a).sum() / s) for a in (X, Y, Z)]


def run(case, T=60.0, lam=1.0, sig=1.5, rho_ref=0.017, E=0.0, gate=True, h=0.5, L=6.0, dt=0.0015, kf=5, Tf=0.6, seed=0,
        knot_R=8, every=1.0):
    os.makedirs(OUT, exist_ok=True)
    u = t10.settled(h, L).copy()
    n = u.shape[1]
    ax = (np.arange(n) - (n - 1) / 2) * h
    if case in ("knot", "knot_grad", "plain_grad"):
        pass
    elif case == "noknot":
        pass
    elif case == "seed":                       # knot + small W seed
        R = np.sqrt(sum(a ** 2 for a in np.meshgrid(ax, ax, ax, indexing="ij")))
        u = g9.init(h, L)[2] * 0 + u          # placeholder to keep shape
        u[2] = 0.5 * np.exp(-(R / 1.0) ** 2)
    lift = s12.hedgehog(n, knot_R) if case in ("knot", "seed", "knot_grad") else np.full((n, n, n), s12.ID_LIFT)
    fd = KnotField(lift, lam, Tf=Tf, seed=seed)
    use_gate = gate and case != "plain_grad"
    Xg, Yg, Zg = np.meshgrid(ax, ax, ax, indexing="ij")
    zb = ax
    bc = [e.GINF + E * Zg[0], e.GINF + E * Zg[-1], e.GINF + E * Zg[:, 0], e.GINF + E * Zg[:, -1],
          np.full(Zg[:, :, 0].shape, e.GINF + E * (zb[0] - h)), np.full(Zg[:, :, -1].shape, e.GINF + E * (zb[-1] + h))]
    if E != 0:                                    # start G from the linear background plus the particle deficit
        u[3] = u[3] + E * Zg
    rows = []
    t0 = time.time()
    nsteps = int(round(T / dt))
    gfac = np.ones((n, n, n))
    for i in range(nsteps + 1):
        if i % kf == 0:
            rho = RHO[fd.lift]
            if use_gate:
                rs = gaussian_filter(rho, sig / h, mode="constant")
                gfac = np.clip(rs / rho_ref, 0, 1)
        if i % int(round(every / dt)) == 0:
            rho = RHO[fd.lift]
            rows.append(dict(t=round(i * dt, 4), M=float(np.sum(K * u[3] * u[2]) * h ** 3), Wmax=float(u[2].max()),
                             Wcen=centre(u[2], ax), Kcen=centre(rho, ax), deg=s12.degree(fd.lift, PTS[0]),
                             E=fd.energy(), gmax=float(gfac.max()), wall=round(time.time() - t0, 1)))
            if not np.isfinite(u).all():
                break
        if i == nsteps:
            break
        if i % kf == 0:
            fd.Wf = u[2]
            fd.sweep()
        X, Y, W, G = u
        s = K * G * W
        u = u + dt * np.array([g9.lap(X, h) + e.A - (1 + e.B) * X + X * X * Y + s,
                               e.DY * g9.lap(Y, h) + e.B * X - X * X * Y,
                               e.DW * g9.lap(W, h) - W + e.C * Y * W * W * gfac - W ** 3,
                               lap_bc(G, h, bc) - s])
    tag = f"{case}_lam{lam:g}_sig{sig:g}_E{E:g}_s{seed}"
    prof = g9.gprofile(np.sqrt(Xg ** 2 + Yg ** 2 + Zg ** 2), u[3] - E * Zg, h)
    rec = dict(case=case, T=T, lam=lam, sig=sig, E=E, gate=use_gate, h=h, L=L, dt=dt, kf=kf, Tf=Tf, seed=seed,
               knot_R=knot_R, rho_ref=rho_ref, rows=rows, Gprofile=prof)
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    r0, r1 = rows[0], rows[-1]
    print("RESULT", tag, "M", round(r0["M"], 3), "->", round(r1["M"], 3), "Wmax", round(r1["Wmax"], 2), "deg",
          r0["deg"], "->", r1["deg"], "Wcen", [round(x, 2) for x in r1["Wcen"]], "Kcen",
          [round(x, 2) for x in r1["Kcen"]], "wall", round(time.time() - t0), flush=True)
    return rec


PTS = [np.array([-0.3, 0.5, -0.6, 0.55]) / np.linalg.norm([-0.3, 0.5, -0.6, 0.55])]

if __name__ == "__main__":
    a = sys.argv
    kw = {}
    for x in a[2:]:
        k, v = x.split("=")
        kw[k] = int(v) if k in ("seed", "kf", "knot_R") else (v == "1" if k == "gate" else float(v))
    run(a[1], **kw)
