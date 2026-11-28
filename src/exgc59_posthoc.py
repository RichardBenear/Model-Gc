"""EX-GC59 post-hoc diagnostics (not scored), after the scored runs were seen:
(1) probe response to uniform-gradient shifts in flat N: rigid rotation beta = s(-y, x), xy shear beta = s(y, x), diagonal shear s(x, -y);
    PS and links only; slope of the odd Sigma angle per unit s;
(2) the in-plane d = 16 and the axis gradient-null runs repeated, tracking the angle of the probe's TOTAL angular momentum
    J = angmom(P) + 1/2 <Sigma> in the xy plane, besides the Sigma angle;
(3) the in-plane Sigma rate predicted from the lattice gradients of beta at the probe with the response coefficients from (1).
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G, exgc53 as T53, exgc57 as E57, exgc58 as X, exgc59 as E


def run_track(Ns, Bs, rule):
    psi = E.probe(); C = X.curl_op(Bs) if Bs is not None else None; out = dict(phi=[], phiJ=[])
    def H(p):
        o = G.Hm(p, Ns, E.M, 0.0)
        if Bs is not None:
            o = o + T53.Hb(p, Bs) + (X.Hpl(p, C) if rule == "PS" else 0)
        return o
    def rec(p):
        sx, sy = E.spin_xy(p); L = T53.angmom(T53.link_P(p)); out["phi"].append(np.arctan2(sy, sx)); out["phiJ"].append(np.arctan2(L[1] + 0.5 * sy, L[0] + 0.5 * sx))
    rec(psi); t = 0.0
    while t < E.T - 1e-9:
        k1 = -1j * H(psi); k2 = -1j * H(psi + E.DT / 2 * k1); k3 = -1j * H(psi + E.DT / 2 * k2); k4 = -1j * H(psi + E.DT * k3)
        psi = psi + E.DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += E.DT; rec(psi)
    return out


def odd(a, b, k):
    return E.slope((np.array(a[k]) - np.array(b[k])) / 2)


if __name__ == "__main__":
    res = {}; LP = T53.link_pos(E.NB); one = np.ones((E.NB,) * 3); s0 = 1e-5
    def uni(kind):
        B = np.zeros((3,) + one.shape); x, y = LP[0][0], LP[0][1]; xb, yb = LP[1][0], LP[1][1]
        if kind == "rot":
            B[0] = -s0 * y; B[1] = s0 * xb
        elif kind == "shear":
            B[0] = s0 * y; B[1] = s0 * xb
        else:
            B[0] = s0 * x; B[1] = -s0 * yb
        return B
    res["uniform"] = {}
    for kind in ("rot", "shear", "diag"):
        for rule in ("PS", "P"):
            a = run_track(one, uni(kind), rule); b = run_track(one, -uni(kind), rule)
            res["uniform"][f"{kind}_{rule}"] = dict(Sigma=odd(a, b, "phi") / s0, J=odd(a, b, "phiJ") / s0); print("uniform", kind, rule, res["uniform"][f"{kind}_{rule}"], flush=True)
    n = 97; src, kap, N, Nd = E.source_terms(n)
    for key in ("plane_16", "axis_16"):
        x0 = E.pos(n, key); Ns = E.cut_site(N, x0); Bu = E.shift(src["up"]["PS"], kap, "zero"); Bd = E.shift(src["down"]["PS"], kap, "zero")
        a = run_track(Ns, E57.sub(Bu, x0), "PS"); b = run_track(Ns, E57.sub(Bd, x0), "PS"); pr, _ = E.omega_pred(src["up"]["PS"], kap, n, x0)
        # lattice gradients of beta at the probe site (link (x, x+e_a) centred at x + e_a/2): d_y beta_x and d_x beta_y
        x, y, z = x0; sa = 0.5 * (Bu[0][x, y + 1, z] - Bu[0][x, y - 1, z] + Bu[0][x - 1, y + 1, z] - Bu[0][x - 1, y - 1, z]) / 2
        sb = 0.5 * (Bu[1][x + 1, y, z] - Bu[1][x - 1, y, z] + Bu[1][x + 1, y - 1, z] - Bu[1][x - 1, y - 1, z]) / 2
        kr = res["uniform"]["rot_PS"]["Sigma"]; ks = res["uniform"]["shear_PS"]["Sigma"]
        row = dict(Omega_pred=pr, Sigma=odd(a, b, "phi"), J=odd(a, b, "phiJ"), dy_bx=float(sa), dx_by=float(sb), Sigma_pred_decomp=float(kr * (sb - sa) / 2 + ks * (sa + sb) / 2))
        row["Sigma_ratio"] = row["Sigma"] / pr; row["J_ratio"] = row["J"] / pr; row["decomp_ratio"] = row["Sigma"] / row["Sigma_pred_decomp"]
        if key == "axis_16":
            Bsub = E57.sub(Bu, x0); amp = float(np.abs(Bsub[:, 12:36, 12:36, 12:36]).max()) * 8.0 * np.exp(0.5); Bg = E57.grad_field(n, x0, amp)
            g = run_track(Ns, E57.sub(Bg, x0), "PS"); g2 = run_track(Ns, -E57.sub(Bg, x0), "PS")
            row["grad"] = dict(Sigma_rel=odd(g, g2, "phi") / row["Sigma"], J_rel=odd(g, g2, "phiJ") / row["J"])
        res[key] = row; print(key, json.dumps(row, default=float), flush=True)
    (E.OUT / "posthoc.json").write_text(json.dumps(res, indent=1, default=float))
