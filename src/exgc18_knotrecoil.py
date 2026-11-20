"""EX-GC18 (exploratory): the EX-GC16/17 recoil rule on the EX-GC14/15 knot particle.

Deterministic Gc1 chemistry plus the knot frame field (EX-GC14 coupling: W autocatalysis gated by smoothed knot
density, knot pulled to W with lam). Recoil (mean-field form of EX-GC16/17): the overdamped recoil velocity is
    v = -2 D_G sum(s grad ln G) / sum(s),   s = K G W,   D_G = 1,
evaluated on the actual field every frame sweep. The displacement accumulates; each time it reaches h along an axis
the KNOT is shifted one site (mode 'knot'; the chemistry must follow through the gate), or all fields are shifted
rigidly (mode 'plain', no knot). Our construction.
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
import exgc10_tether as t10  # noqa: E402
import exgc12_soliton as s12  # noqa: E402
import exgc14_knotgc as k14  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC18")


def settled(K, h, L):
    if K == 0.2:
        return t10.settled(h, L)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"settled_K{K:g}_h{h:g}_L{L:g}.npy")
    if os.path.exists(path):
        return np.load(path)
    g9.K = K
    ax, R, u, rows = g9.deterministic(h, L, 30.0, 0.0015 * min(1.0, (h / 0.5) ** 2))
    g9.K = 0.2
    print("settled K", K, "M", rows[-1]["M"], "W0", rows[-1]["W0"], flush=True)
    np.save(path, u)
    return u


def shift(a, k, axis, fill):
    out = np.roll(a, k, axis=axis)
    sl = [slice(None)] * a.ndim
    sl[axis] = slice(0, k) if k > 0 else slice(k, None)
    out[tuple(sl)] = fill
    return out


def run(mode, K=0.2, E=0.02, T=30.0, recoil=True, lam=1.0, sig=1.5, rho_ref=0.017, h=0.25, L=6.0, kf=20, Tf=0.6,
        seed=0, knot_R=8, every=1.0, field=False):
    """field=True (EX-GC18b): recoil as local advection of W, dW/dt += -div(W v), v = -2 D_G grad ln G."""
    os.makedirs(OUT, exist_ok=True)
    dt = 0.0015 * min(1.0, (h / 0.5) ** 2)
    u = settled(K, h, L).copy()
    n = u.shape[1]
    ax = (np.arange(n) - (n - 1) / 2) * h
    Xg, Yg, Zg = np.meshgrid(ax, ax, ax, indexing="ij")
    Gbg = e.GINF + E * Zg
    u[3] = u[3] + E * Zg
    bc = [e.GINF + E * Zg[0], e.GINF + E * Zg[-1], e.GINF + E * Zg[:, 0], e.GINF + E * Zg[:, -1],
          np.full(Zg[:, :, 0].shape, e.GINF + E * (ax[0] - h)), np.full(Zg[:, :, -1].shape, e.GINF + E * (ax[-1] + h))]
    knot = mode == "knot"
    lift = s12.hedgehog(n, knot_R) if knot else np.full((n, n, n), s12.ID_LIFT)
    fd = k14.KnotField(lift, lam, Tf=Tf, seed=seed)
    gfac = np.ones((n, n, n))
    disp = np.zeros(3)
    shifts = np.zeros(3, int)
    vlast = np.zeros(3)
    rows = []
    t0 = time.time()
    nsteps = int(round(T / dt))
    for i in range(nsteps + 1):
        if i % kf == 0 and knot:
            rs = gaussian_filter(k14.RHO[fd.lift], sig / h, mode="constant")
            gfac = np.clip(rs / rho_ref, 0, 1)
        if i % int(round(every / dt)) == 0:
            rho = k14.RHO[fd.lift]
            rows.append(dict(t=round(i * dt, 4), M=float(np.sum(K * u[3] * u[2]) * h ** 3), Wmax=float(u[2].max()),
                             Wcen=k14.centre(u[2], ax), Kcen=k14.centre(rho, ax) if knot else None,
                             deg=s12.degree(fd.lift, k14.PTS[0]) if knot else None, v=vlast.tolist(),
                             shifts=shifts.tolist(), gmax=float(gfac.max()), wall=round(time.time() - t0, 1)))
            if not np.isfinite(u).all() or u[2].max() < 1e-3:
                break
        if i == nsteps:
            break
        if i % kf == 0:
            if knot:
                fd.Wf = u[2]
                fd.sweep()
            if recoil and not field:
                G = u[3]
                s = K * G * u[2]
                lg = np.log(np.maximum(G, 1e-9))
                gr = np.gradient(lg, h)
                ss = s.sum()
                vlast = np.array([-2.0 * (s * g).sum() / ss for g in gr]) if ss > 0 else np.zeros(3)
                disp += vlast * dt * kf
                for a in range(3):
                    while abs(disp[a]) >= h:
                        k = int(np.sign(disp[a]))
                        disp[a] -= k * h
                        shifts[a] += k
                        if knot:
                            fd.lift = shift(fd.lift, k, a, s12.ID_LIFT)
                        else:
                            u[0] = shift(u[0], k, a, e.A)
                            u[1] = shift(u[1], k, a, e.B / e.A)
                            u[2] = shift(u[2], k, a, 0.0)
                            u[3] = Gbg + shift(u[3] - Gbg, k, a, 0.0)
        X, Y, W, G = u
        s = K * G * W
        adv = 0.0
        if recoil and field:
            vf = [-2.0 * g for g in np.gradient(np.log(np.maximum(G, 1e-9)), h)]
            adv = -sum(np.gradient(W * vf[a], h, axis=a) for a in range(3))
            if i % kf == 0:
                ws = W.sum()
                vlast = np.array([(W * v).sum() / ws for v in vf]) if ws > 0 else np.zeros(3)
        u = u + dt * np.array([g9.lap(X, h) + e.A - (1 + e.B) * X + X * X * Y + s,
                               e.DY * g9.lap(Y, h) + e.B * X - X * X * Y,
                               e.DW * g9.lap(W, h) - W + e.C * Y * W * W * gfac - W ** 3 + adv,
                               k14.lap_bc(G, h, bc) - s])
    tt = np.array([r["t"] for r in rows])
    wz = np.array([r["Wcen"][2] for r in rows])
    sel = tt >= 10
    drift = float(np.polyfit(tt[sel], wz[sel], 1)[0]) if sel.sum() > 3 else float("nan")
    tag = f"{mode}_K{K:g}_E{E:g}_rec{int(recoil)}_s{seed}" + ("_field" if field else "")
    rec = dict(mode=mode, K=K, E=E, T=T, recoil=recoil, lam=lam, sig=sig, rho_ref=rho_ref, h=h, L=L, dt=dt, kf=kf, field=field,
               Tf=Tf, seed=seed, rows=rows, drift_z=drift, platform=sys.version.split()[0] + " numpy " + np.__version__)
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    r0, r1 = rows[0], rows[-1]
    print("RESULT", tag, "t", r1["t"], "M", round(r0["M"], 3), "->", round(r1["M"], 3), "drift_z", round(drift, 5),
          "v", np.round(r1["v"], 4).tolist(), "shifts", r1["shifts"], "Wcen", np.round(r1["Wcen"], 2).tolist(),
          "Kcen", (np.round(r1["Kcen"], 2).tolist() if knot else None), "deg", r1["deg"], "wall",
          round(time.time() - t0), flush=True)
    return rec


if __name__ == "__main__":
    a = sys.argv
    kw = {}
    for x in a[2:]:
        k, v = x.split("=")
        kw[k] = int(v) if k in ("seed", "kf", "knot_R") else (v == "1" if k in ("recoil", "field") else float(v))
    run(a[1], **kw)
