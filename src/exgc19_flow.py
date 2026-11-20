"""EX-GC19 (exploratory): a flow (momentum) field for the Gc particle.

Deterministic Gc1 chemistry (h = 0.25) plus a velocity field V (3 components):
  * X, Y, W deviations from background are carried by the flow: df/dt += -div((f - f0) V); G is not carried.
  * Inside the core (weight s/(s + eps), s = K G W, mass density proportional to s, coasting time tau):
        dV/dt = [-2 D_G grad ln G - drag * V] / tau            (recoil push per absorption, EX-GC16-18)
  * Cohesion: dV/dt += -kappa (V - <V>), <V> = s-weighted Gaussian average (width sig_v).
Optional knot (EX-GC14 coupling). Our construction.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc09_lattice as g9  # noqa: E402
import exgc12_soliton as s12  # noqa: E402
import exgc14_knotgc as k14  # noqa: E402
import exgc18_knotrecoil as r18  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC19")


def div(F, h):
    return sum(np.gradient(F[a], h, axis=a) for a in range(3))


def run(case, K=0.2, E=0.02, T=30.0, drag=1.0, kappa=1.0, tau=20.0, kick=0.0, knot=False, sig_v=1.5, h=0.25, L=6.0,
        kf=20, lam=1.0, sig=1.5, rho_ref=0.017, Tf=0.6, seed=0, knot_R=8, every=1.0):
    os.makedirs(OUT, exist_ok=True)
    dt = 0.0015 * min(1.0, (h / 0.5) ** 2)
    u = r18.settled(K, h, L).copy()
    n = u.shape[1]
    ax = (np.arange(n) - (n - 1) / 2) * h
    Xg, Yg, Zg = np.meshgrid(ax, ax, ax, indexing="ij")
    u[3] = u[3] + E * Zg
    bc = [e.GINF + E * Zg[0], e.GINF + E * Zg[-1], e.GINF + E * Zg[:, 0], e.GINF + E * Zg[:, -1],
          np.full(Zg[:, :, 0].shape, e.GINF + E * (ax[0] - h)), np.full(Zg[:, :, -1].shape, e.GINF + E * (ax[-1] + h))]
    f0 = np.array([e.A, e.B / e.A, 0.0])
    V = np.zeros((3, n, n, n))
    V[2] = -kick
    lift = s12.hedgehog(n, knot_R) if knot else np.full((n, n, n), s12.ID_LIFT)
    fd = k14.KnotField(lift, lam, Tf=Tf, seed=seed)
    gfac = np.ones((n, n, n))
    rows = []
    t0 = time.time()
    nsteps = int(round(T / dt))
    for i in range(nsteps + 1):
        X, Y, W, G = u
        s = K * G * W
        if i % kf == 0:
            if knot:
                rs = gaussian_filter(k14.RHO[fd.lift], sig / h, mode="constant")
                gfac = np.clip(rs / rho_ref, 0, 1)
            wgt = s / (s + 1e-3 * s.max())
            push = [-2.0 * g for g in np.gradient(np.log(np.maximum(G, 1e-9)), h)]
            sm = gaussian_filter(s, sig_v / h, mode="constant") + 1e-12
            Vbar = np.array([gaussian_filter(s * V[a], sig_v / h, mode="constant") / sm for a in range(3)])
            dV = np.array([wgt * (push[a] - drag * V[a]) / tau for a in range(3)]) - kappa * (V - Vbar)
            V = V + kf * dt * dV
        if i % int(round(every / dt)) == 0:
            ss = s.sum()
            Vc = [float((s * V[a]).sum() / ss) for a in range(3)] if ss > 0 else [0.0] * 3
            sp = float(np.sqrt(((s * ((V - np.array(Vc)[:, None, None, None]) ** 2).sum(0)).sum()) / ss)) if ss > 0 else 0
            rows.append(dict(t=round(i * dt, 4), M=float(s.sum() * h ** 3), Wmax=float(W.max()), Wcen=k14.centre(W, ax),
                             Kcen=k14.centre(k14.RHO[fd.lift], ax) if knot else None,
                             deg=s12.degree(fd.lift, k14.PTS[0]) if knot else None, Vcore=Vc, Vspread=sp,
                             Vmax=float(np.abs(V).max()), wall=round(time.time() - t0, 1)))
            if not np.isfinite(u).all() or W.max() < 1e-3 or not np.isfinite(V).all():
                break
        if i == nsteps:
            break
        if knot and i % kf == 0:
            fd.Wf = W
            fd.sweep()
        adv = [-div((u[j] - f0[j]) * V, h) for j in range(3)]
        u = u + dt * np.array([g9.lap(X, h) + e.A - (1 + e.B) * X + X * X * Y + s + adv[0],
                               e.DY * g9.lap(Y, h) + e.B * X - X * X * Y + adv[1],
                               e.DW * g9.lap(W, h) - W + e.C * Y * W * W * gfac - W ** 3 + adv[2],
                               k14.lap_bc(G, h, bc) - s])
    rec = dict(case=case, K=K, E=E, T=T, drag=drag, kappa=kappa, tau=tau, kick=kick, knot=knot, sig_v=sig_v, h=h, L=L,
               dt=dt, kf=kf, lam=lam, rows=rows, platform=sys.version.split()[0] + " numpy " + np.__version__)
    json.dump(rec, open(os.path.join(OUT, case + ".json"), "w"), indent=1)
    r0, r1 = rows[0], rows[-1]
    print("RESULT", case, "t", r1["t"], "M", round(r0["M"], 3), "->", round(r1["M"], 3), "Wcen_z", round(r1["Wcen"][2], 3),
          "Vcore_z", round(r1["Vcore"][2], 4), "Vspread", round(r1["Vspread"], 4), "Vmax", round(r1["Vmax"], 3),
          "Kcen_z", (round(r1["Kcen"][2], 3) if knot else None), "wall", round(time.time() - t0), flush=True)
    return rec


if __name__ == "__main__":
    a = sys.argv
    kw = {}
    for x in a[2:]:
        k, v = x.split("=")
        kw[k] = int(v) if k in ("seed", "kf", "knot_R") else (v == "1" if k == "knot" else float(v))
    run(a[1], **kw)
