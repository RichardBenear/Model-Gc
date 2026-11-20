"""EX-GC16 summary: compare measured recoil forces with the mean-field lattice prediction.

Separate path from the simulation: the mean token field is the steady solution of the lattice diffusion equation
(D = 1, h = 1) with Dirichlet faces G = 1 + E z and a sink of m0 tokens per time at each consumer site. The predicted
recoil on a unit is m0 * (net direction of arrivals), with arrivals in direction e proportional to the token count of
the neighbour behind it. Shares only the geometry (make_units) with the simulator.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc16_recoil as r  # noqa: E402  (geometry only)


def mean_field(units, L, Om, E):
    N = 2 * L + 1
    idx = np.arange(N ** 3).reshape(N, N, N)
    z = np.arange(N) - L
    rows, cols, vals = [], [], []
    b = np.zeros(N ** 3)
    for a in range(3):
        for s in (1, -1):
            sl_in = [slice(None)] * 3
            sl_out = [slice(None)] * 3
            if s == 1:
                sl_in[a], sl_out[a] = slice(0, -1), slice(1, None)
            else:
                sl_in[a], sl_out[a] = slice(1, None), slice(0, -1)
            rows.append(idx[tuple(sl_in)].ravel())
            cols.append(idx[tuple(sl_out)].ravel())
            vals.append(np.ones(rows[-1].size))
            face = [slice(None)] * 3
            face[a] = N - 1 if s == 1 else 0
            Z = np.broadcast_to(z[None, None, :], (N, N, N))[tuple(face)]
            if a == 2:
                Zb = np.full_like(Z, (L + 1) if s == 1 else (-L - 1), dtype=float)
            else:
                Zb = Z.astype(float)
            b[idx[tuple(face)].ravel()] -= Om * (1 + E * Zb).ravel()
    rows.append(idx.ravel())
    cols.append(idx.ravel())
    vals.append(-6.0 * np.ones(N ** 3))
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N ** 3, N ** 3))
    for u in units:
        for s in u["sites"]:
            b[idx[tuple(s)]] += u["rate"] / len(u["sites"])          # spread units: sink shared (approximation)
    return spla.spsolve(A.tocsc(), b).reshape(N, N, N)


def pred_force(units, n):
    out = []
    for u in units:
        F = np.zeros(3)
        tot = 0.0
        for s in u["sites"]:
            for k, (a, sg) in enumerate(r.DIRS):
                src = np.array(s)
                src[a] -= sg                                         # token travelling sg along a came from here
                F += n[tuple(src)] * r.EVEC[k]
                tot += n[tuple(src)]
        out.append(u["rate"] * F / tot)
    return np.array(out)


def main():
    files = sorted(glob.glob(os.path.join(r.OUT, "*_s*.json")))
    groups = defaultdict(list)
    for f in files:
        rec = json.load(open(f))
        groups[(rec["kind"], rec["m0"], rec["E"], rec["d"], rec["L"])].append(rec)
    rows = []
    for (kind, m0, E, d, L), recs in sorted(groups.items()):
        units, B = r.make_units(kind, L, m0, d)
        n = mean_field(units, L, recs[0]["Om"], E)
        Fp = pred_force(units, n)
        fe = np.array([x["fall_exp"] for x in recs])
        fu = np.array([x["force_exp"] for x in recs])                 # seeds x units x 3
        T = sum(x["T"] - x["burn"] for x in recs)
        row = dict(kind=kind, m0=m0, E=E, d=d, L=L, seeds=len(recs), T_total=T,
                   fall_z=float(fe[:, 2].mean()), fall_z_sem=float(fe[:, 2].std(ddof=1) / np.sqrt(len(recs)))
                   if len(recs) > 1 else float("nan"),
                   fall_xy_rms=float(np.sqrt((fe[:, :2] ** 2).mean())),
                   pred_fall=(Fp.sum(0) / B).tolist(), ideal_fall_z=-m0 * E / 3,
                   unit_force=fu.mean(0).tolist(), pred_unit_force=Fp.tolist(),
                   absorbed=np.mean([x["absorbed_rate"] for x in recs], 0).tolist())
        rows.append(row)
        print(f"{kind:12s} m0={m0:<5g} E={E:<5g} d={d} L={L} seeds={len(recs)}  fall_z={row['fall_z']:+.4f}"
              f" ±{row['fall_z_sem']:.4f}  pred={row['pred_fall'][2]:+.4f}  ideal={row['ideal_fall_z']:+.4f}"
              f"  xy_rms={row['fall_xy_rms']:.4f}")
        if E == 0 and len(units) == 2:
            a = fu[:, 0, :] - fu[:, 1, :]
            ax = 0 if kind == "pair_x" else 2
            meas = a[:, ax] / 2
            print(f"    pair pull (unit0 toward unit1, along axis {ax}): {meas.mean():+.4f}"
                  f" ±{(meas.std(ddof=1) / np.sqrt(len(meas)) if len(meas) > 1 else float('nan')):.4f}"
                  f"   pred {(Fp[0, ax] - Fp[1, ax]) / 2:+.4f}")
            row["pair_pull"] = float(meas.mean())
            row["pair_pull_pred"] = float((Fp[0, ax] - Fp[1, ax]) / 2)
    json.dump(rows, open(os.path.join(r.OUT, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
