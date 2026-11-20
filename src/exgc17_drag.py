"""EX-GC17 (exploratory): drag and acceleration response of EX-GC16 recoil consumers.

Same token lattice and consumer rule as EX-GC16 (E = 0). The particle is held fixed while the token gas drifts past it
with velocity u(t) along z: per-step hop probabilities dt along x, y and dt (1 +- u/2) along +-z. In the particle's
frame this is the particle moving at v = -u through the gas. Modes: steady u, or u = u0 sin(w t), for which the force
is projected on sin and cos to give the in-phase (drag) and quadrature (acceleration) parts. Our construction.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc16_recoil as r16  # noqa: E402

ROOT = r16.ROOT
OUT = os.path.join(ROOT, "results", "EX-GC17")
DIRS, EVEC = r16.DIRS, r16.EVEC


def probs(dt, u):
    return np.array([dt, dt, dt, dt, dt * (1 + u / 2), dt * (1 - u / 2)])


def hop(n, p, rng, Om):
    N = n.shape[0]
    leave = rng.binomial(n, p.sum())
    new = n - leave
    rem = leave
    left = p.sum()
    shares = []
    for k in range(5):
        x = rng.binomial(rem, min(1.0, p[k] / left))
        shares.append(x)
        rem = rem - x
        left -= p[k]
    shares.append(rem)
    arrivals = []
    for k, (axis, sgn) in enumerate(DIRS):
        s = shares[k]
        moved = np.zeros_like(n)
        src = [slice(None)] * 3
        dst = [slice(None)] * 3
        if sgn == 1:
            src[axis], dst[axis] = slice(0, -1), slice(1, None)
        else:
            src[axis], dst[axis] = slice(1, None), slice(0, -1)
        moved[tuple(dst)] = s[tuple(src)]
        face = [slice(None)] * 3
        face[axis] = 0 if sgn == 1 else N - 1
        moved[tuple(face)] += rng.poisson(Om * p[k], size=(N, N))
        new += moved
        arrivals.append(moved)
    return new, arrivals


def run(kind, m0=500.0, u=0.1, w=0.0, d=4, L=10, Om=1000.0, dt=0.1, T=1000.0, burn=50.0, seed=0):
    """w = 0: steady drift u. w > 0: u(t) = u sin(w t)."""
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(seed + 1000)
    N = 2 * L + 1
    n = rng.poisson(Om, size=(N, N, N)).astype(np.int64)
    units, B = r16.make_units(kind, L, m0, d)
    nsteps = int(round(T / dt))
    nburn = int(round(burn / dt))
    F = np.zeros(3)
    Fs = np.zeros(3)
    Fc = np.zeros(3)
    S2 = C2 = SC = 0.0
    absorbed = 0
    t0 = time.time()
    for i in range(nsteps):
        t = (i + 0.5) * dt
        ut = u * np.sin(w * t) if w > 0 else u
        n, arr = hop(n, probs(dt, ut), rng, Om)
        Fstep = np.zeros(3)
        for un in units:
            A = np.array([[arr[k][tuple(s)] for k in range(6)] for s in un["sites"]])
            tot = int(A.sum())
            kabs = min(int(rng.poisson(un["rate"] * dt)), tot)
            if kabs == 0:
                continue
            take = rng.multivariate_hypergeometric(A.ravel(), kabs).reshape(A.shape)
            for j, s in enumerate(un["sites"]):
                n[tuple(s)] -= int(take[j].sum())
            Fstep += kabs * (A.sum(0) @ EVEC) / tot
            if i >= nburn:
                absorbed += kabs
        if i >= nburn:
            F += Fstep
            if w > 0:
                sn, cs = np.sin(w * t), np.cos(w * t)
                Fs += Fstep * sn
                Fc += Fstep * cs
                S2 += sn * sn
                C2 += cs * cs
                SC += sn * cs
    Tm = (nsteps - nburn) * dt
    rec = dict(kind=kind, m0=m0, u=u, w=w, d=d, L=L, Om=Om, dt=dt, T=T, burn=burn, seed=seed, B=B,
               force=(F / Tm).tolist(), absorbed_rate=absorbed / Tm, wall=time.time() - t0)
    if w > 0:
        # least squares F = a sin + b cos (+ mean, small): solve the 2x2 normal equations per component
        M = np.array([[S2, SC], [SC, C2]])
        ab = np.linalg.solve(M, np.vstack([Fs, Fc]))          # 2 x 3 (per unit time-step sums -> per step force)
        a, b = ab[0] / dt, ab[1] / dt                           # force per unit time
        # u = u0 sin; particle velocity v = -u, acceleration acc = -u0 w cos; F = -gamma v + kappa acc
        rec["gamma_per_B"] = float(a[2] / u / B)
        rec["kappa_per_B"] = float(-b[2] / (u * w) / B)
        rec["phase_fit"] = dict(sin=a.tolist(), cos=b.tolist())
    else:
        rec["gamma_per_B"] = float(F[2] / Tm / u / B)
    tag = f"{kind}_m{m0:g}_u{u:g}_w{w:g}_s{seed}"
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    print("RESULT", tag, "force", np.round(rec["force"], 4).tolist(), "gamma/B", round(rec["gamma_per_B"], 3),
          "(weak-limit m0/6 =", round(m0 / 6, 3), ")", "kappa/B", round(rec.get("kappa_per_B", float("nan")), 3),
          "abs", round(rec["absorbed_rate"], 1), "wall", round(rec["wall"]), flush=True)
    return rec


if __name__ == "__main__":
    a = sys.argv
    kw = {}
    for x in a[2:]:
        k, v = x.split("=")
        kw[k] = int(v) if k in ("d", "L", "seed") else float(v)
    run(a[1], **kw)
