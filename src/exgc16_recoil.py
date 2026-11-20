"""EX-GC16 (exploratory): recoil rule for universal free fall, on a G token lattice.

G tokens random-walk on a cubic lattice (h = 1) with a reservoir at the faces holding G = 1 + E z. A particle is made
of B consumer units; each unit absorbs tokens at a FIXED rate m0 (zero order in G), drawn from the tokens arriving at
its site(s) in that step, and each absorbed token pushes the particle by one unit along the token's direction of travel.
Particles are held fixed; we record the push per unit time (force). 'Fall rate' = force / B. Our construction.
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC16")
DIRS = ((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1))      # (axis, sign) of travel
EVEC = np.array([[s if a == ax else 0 for a in range(3)] for ax, s in DIRS], float)


def hop(n, p_hop, rng, L, Om, E):
    """Tokens leave with prob p_hop (6 directions). Returns new grid and per-direction moved arrays (arrivals)."""
    N = n.shape[0]
    leave = rng.binomial(n, p_hop)
    new = n - leave
    rem = leave.copy()
    shares = []
    for k in range(5):
        x = rng.binomial(rem, 1.0 / (6 - k))
        shares.append(x)
        rem = rem - x
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
        # inflow from the reservoir through the face where tokens travel inward in this direction
        face = [slice(None)] * 3
        face[axis] = 0 if sgn == 1 else N - 1
        z = np.arange(N) - L
        if axis == 2:
            zout = (-L - 1) if sgn == 1 else (L + 1)
            Gout = np.full((N, N), 1.0 + E * zout)
        else:
            Gout = np.broadcast_to(1.0 + E * z[None, :], (N, N))
        moved[tuple(face)] += rng.poisson(Om * Gout * p_hop / 6.0)
        new += moved
        arrivals.append(moved)
    return new, arrivals


def make_units(kind, L, m0, d=4):
    c = np.array([L, L, L])
    nb7 = [c + o for o in ([0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    if kind == "single":
        return [dict(sites=[c], rate=m0)], 1
    if kind == "double_site":
        return [dict(sites=[c], rate=2 * m0)], 2
    if kind == "spread7":
        return [dict(sites=nb7, rate=m0)], 1
    if kind == "pair_z":
        return [dict(sites=[c + [0, 0, -d // 2]], rate=m0), dict(sites=[c + [0, 0, d - d // 2]], rate=m0)], 2
    if kind == "pair_x":
        return [dict(sites=[c + [-(d // 2), 0, 0]], rate=m0), dict(sites=[c + [d - d // 2, 0, 0]], rate=m0)], 2
    raise ValueError(kind)


def run(kind, m0=50.0, E=0.02, d=4, L=10, Om=1000.0, dt=0.1, T=200.0, burn=50.0, seed=0):
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(seed)
    N = 2 * L + 1
    z = np.arange(N) - L
    n = rng.poisson(Om * (1.0 + E * z)[None, None, :] * np.ones((N, N, N))).astype(np.int64)
    units, B = make_units(kind, L, m0, d)
    p_hop = 6 * dt                                    # D = 1, h = 1: per-direction hop prob dt
    F_real = np.zeros((len(units), 3))
    F_exp = np.zeros((len(units), 3))
    absorbed = np.zeros(len(units))
    nsteps = int(round(T / dt))
    nburn = int(round(burn / dt))
    t0 = time.time()
    Gsite = []
    for i in range(nsteps):
        n, arr = hop(n, p_hop, rng, L, Om, E)
        for u_i, u in enumerate(units):
            A = np.array([[arr[k][tuple(s)] for k in range(6)] for s in u["sites"]])      # sites x dirs
            tot = int(A.sum())
            kabs = min(int(rng.poisson(u["rate"] * dt)), tot)
            if kabs == 0:
                continue
            take = rng.multivariate_hypergeometric(A.ravel(), kabs).reshape(A.shape)
            for j, s in enumerate(u["sites"]):
                n[tuple(s)] -= int(take[j].sum())
            if i >= nburn:
                F_real[u_i] += take.sum(0) @ EVEC
                F_exp[u_i] += kabs * (A.sum(0) @ EVEC) / tot
                absorbed[u_i] += kabs
        if i >= nburn and i % 10 == 0:
            Gsite.append([float(n[tuple(u["sites"][0])]) / Om for u in units])
    Tm = (nsteps - nburn) * dt
    rec = dict(kind=kind, m0=m0, E=E, d=d, L=L, Om=Om, dt=dt, T=T, burn=burn, seed=seed, B=B,
               force_real=(F_real / Tm).tolist(), force_exp=(F_exp / Tm).tolist(),
               absorbed_rate=(absorbed / Tm).tolist(), G_at_unit=np.mean(Gsite, axis=0).tolist(),
               fall_exp=((F_exp.sum(0) / Tm) / B).tolist(), fall_real=((F_real.sum(0) / Tm) / B).tolist(),
               wall=time.time() - t0)
    tag = f"{kind}_m{m0:g}_E{E:g}_d{d}_s{seed}"
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    print("RESULT", tag, "fall_exp", np.round(rec["fall_exp"], 4), "fall_real", np.round(rec["fall_real"], 3),
          "per-unit force_exp", np.round(rec["force_exp"], 4).tolist(), "abs", np.round(rec["absorbed_rate"], 1).tolist(),
          "G_unit", np.round(rec["G_at_unit"], 3).tolist(), "wall", round(rec["wall"]), flush=True)
    return rec


if __name__ == "__main__":
    a = sys.argv
    kw = {}
    for x in a[2:]:
        k, v = x.split("=")
        kw[k] = int(v) if k in ("d", "L", "seed") else float(v)
    run(a[1], **kw)
