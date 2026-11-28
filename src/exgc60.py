"""EX-GC60: what the spin detector responds to, with the coupling unchanged. 60A uniform-gradient fields (rotation, xy shear, sum) at
sigma = 3, 4, 5 in two orientations, tracking the angles of Sigma, L and J = L + Sigma/2; 60B gauge test with the lattice resync coupling
H_g = -(i/2)[chi, H^2]; 60C frozen-calibration prediction at a new diagonal site in 59's source field; 60D 59's in-plane case at dt/2.
Plan: spec/EX-GC60_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc57 as E57
import exgc58 as X
import exgc59 as E

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC60"
M = 0.6; NB = 48; T = 10.0; S0 = 1e-5
CHI = {"x": (1 / np.sqrt(2), 1 / np.sqrt(2), 0, 0), "y": (1 / np.sqrt(2), 1j / np.sqrt(2), 0, 0)}
PLANE = {"xy": (0, 1, "x"), "yz": (1, 2, "y")}                                   # (first, second axis, initial spin direction)


def probe(sig, spin):
    c = NB // 2; psi, _, _ = S.packet_k(NB, M, 0.0, (c, c, c), sig, (0, 0, 0), chi=CHI[spin]); return psi


def angles(psi, plane):
    a, b, _ = PLANE[plane]; s = [float(np.real(np.vdot(psi, psi @ X.SIG[k].T))) for k in range(3)]; L = T53.angmom(T53.link_P(psi))
    J = [L[k] + 0.5 * s[k] for k in range(3)]
    return np.arctan2(s[b], s[a]), np.arctan2(L[b], L[a]), np.arctan2(J[b], J[a]), s


def run(Ns, Bs, rule, sig, plane, dt=0.1, extra=None, psi0=None):
    a, b, sp = PLANE[plane]; psi = probe(sig, sp) if psi0 is None else psi0.copy(); C = X.curl_op(Bs) if (Bs is not None and rule == "PS") else None
    def H(p):
        o = G.Hm(p, Ns, M, 0.0)
        if Bs is not None:
            o = o + T53.Hb(p, Bs) + (X.Hpl(p, C) if C is not None else 0)
        if extra is not None:
            o = o + extra(p)
        return o
    out = dict(S=[], L=[], J=[]); t = 0.0
    def rec(p):
        s, l, j, _ = angles(p, plane); out["S"].append(s); out["L"].append(l); out["J"].append(j)
    rec(psi)
    while t < T - 1e-9:
        k1 = -1j * H(psi); k2 = -1j * H(psi + dt / 2 * k1); k3 = -1j * H(psi + dt / 2 * k2); k4 = -1j * H(psi + dt * k3)
        psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; rec(psi)
    out["dt"] = dt; return out


def slope(phi, dt):
    t = dt * np.arange(len(phi)); A = np.stack([np.ones_like(t), t], 1); return float(np.linalg.lstsq(A, np.asarray(phi), rcond=None)[0][1])


def rates(r, r0, s):
    return {k: slope((np.array(r[k]) - np.array(r0[k])), r["dt"]) / s for k in ("S", "L", "J")}


def field(kind, plane, s=S0):
    a, b, _ = PLANE[plane]; LP = T53.link_pos(NB); B = np.zeros((3, NB, NB, NB))
    xa_on_b = LP[b][a]; xb_on_a = LP[a][b]                                       # coordinate a at b-links, coordinate b at a-links
    if kind in ("rot", "sum"):
        B[a] += -s * xb_on_a; B[b] += s * xa_on_b
    if kind in ("shear", "sum"):
        B[a] += s * xb_on_a; B[b] += s * xa_on_b
    return B


def f_orb(sig):
    p = probe(sig, "x"); c = NB // 2; up, _, _ = S.packet_k(NB, M, 0.0, (c, c, c), sig, (0, 0, 0), chi=T53.UP)
    return 1 - float(np.real(np.vdot(up, up @ X.SIG[2].T)))


# ------------------------------------------------------------------ 60A
def part_A(sig):
    one = np.ones((NB,) * 3); res = dict(sigma=sig, f=f_orb(sig), rows={})
    for plane in ("xy", "yz"):
        r0 = run(one, None, "PS", sig, plane); res["rows"][f"{plane}_zero_again"] = None
        if plane == "xy":
            r0b = run(one, None, "PS", sig, plane); res["C1"] = float(max(np.abs(np.array(r0["S"]) - np.array(r0b["S"]))))
        _, _, _, s0 = angles(probe(sig, PLANE[plane][2]), plane); res[f"S0_{plane}"] = s0
        for kind in ("rot", "shear", "sum"):
            res["rows"][f"{plane}_{kind}_PS"] = rates(run(one, field(kind, plane), "PS", sig, plane), r0, S0)
        if plane == "xy":
            res["rows"]["xy_shear_P"] = rates(run(one, field("shear", plane), "P", sig, plane), r0, S0)
            if sig == 3.0:
                res["rows"]["xy_rot_P"] = rates(run(one, field("rot", plane), "P", sig, plane), r0, S0)
                res["rows"]["xy_shear2_PS"] = rates(run(one, field("shear", plane, 2 * S0), "PS", sig, plane), r0, S0)
        del res["rows"][f"{plane}_zero_again"]
        print("60A sigma", sig, plane, json.dumps({k: v for k, v in res["rows"].items() if k.startswith(plane)}), flush=True)
    return res


# ------------------------------------------------------------------ 60B
def part_B(sig=3.0):
    one = np.ones((NB,) * 3); x = np.arange(NB) - NB // 2; Xg, Yg = np.meshgrid(x, x, indexing="ij"); chi = (-S0 * Xg * Yg)[..., None][:, :, None, :] * np.ones((1, 1, NB, 1))
    Hm = lambda p: G.Hm(p, one, M, 0.0); Hg = lambda p: -0.5j * (chi * Hm(Hm(p)) - Hm(Hm(chi * p)))
    Gop = lambda p: 0.5 * (chi * Hm(p) + Hm(chi * p))
    def expG(p, sgn):
        out = p.copy(); term = p.copy()
        for k in range(1, 12):
            term = (sgn * -1j) / k * Gop(term); out = out + term
        return out
    r0 = run(one, None, "PS", sig, "xy"); rg = run(one, None, "PS", sig, "xy", extra=Hg)
    # gauge-transformed free evolution: psi_free(0) = U^+ psi0, read U psi_free(t)
    psi = expG(probe(sig, "x"), -1); dt = 0.1; t = 0.0; out = dict(S=[], L=[], J=[], dt=dt)
    def rec(p):
        q = expG(p, +1); s, l, j, _ = angles(q, "xy"); out["S"].append(s); out["L"].append(l); out["J"].append(j)
    rec(psi)
    while t < T - 1e-9:
        k1 = -1j * Hm(psi); k2 = -1j * Hm(psi + dt / 2 * k1); k3 = -1j * Hm(psi + dt / 2 * k2); k4 = -1j * Hm(psi + dt * k3); psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; rec(psi)
    rp = run(one, field("shear", "xy"), "PS", sig, "xy")
    res = dict(rate_Hg=rates(rg, r0, S0), rate_gauge=rates(out, r0, S0), rate_P=rates(rp, r0, S0))
    res["Hg_vs_gauge_end"] = float(abs(rg["S"][-1] - out["S"][-1])); res["P_signal_end"] = float(abs(rp["S"][-1] - r0["S"][-1]))
    res["Hg_vs_gauge_rel"] = res["Hg_vs_gauge_end"] / res["P_signal_end"]; res["Hg_over_P"] = res["rate_Hg"]["S"] / res["rate_P"]["S"]
    print("60B", json.dumps(res), flush=True); return res


# ------------------------------------------------------------------ 60C / 60D (59's source field)
def weighted_grads(B, x0, sig=3.0):
    Bs = E57.sub(B, x0); w = np.sum(np.abs(probe(sig, "x")) ** 2, -1); w /= w.sum()
    dyb_x = 0.25 * (np.roll(Bs[0], -1, 1) - np.roll(Bs[0], 1, 1) + np.roll(np.roll(Bs[0], -1, 1), 1, 0) - np.roll(np.roll(Bs[0], 1, 1), 1, 0))
    dxb_y = 0.25 * (np.roll(Bs[1], -1, 0) - np.roll(Bs[1], 1, 0) + np.roll(np.roll(Bs[1], -1, 0), 1, 1) - np.roll(np.roll(Bs[1], 1, 0), 1, 1))
    return float(np.sum(w * dyb_x)), float(np.sum(w * dxb_y))


def part_C(coef):
    n = 97; src, kap, N, Nd = E.source_terms(n); Bu = E.shift(src["up"]["PS"], kap, "zero"); Bd = E.shift(src["down"]["PS"], kap, "zero"); c = n // 2; res = {}
    for key, x0 in (("diag", (c + 12, c, c + 12)), ("axis_16", (c, c, c + 16)), ("plane_16", (c + 16, c, c))):
        sa, sb = weighted_grads(Bu, x0); om = (sb - sa) / 2; sh = (sb + sa) / 2
        pred = {k: coef["rot"][k] * om + coef["shear"][k] * sh for k in ("S", "J")}; row = dict(dy_bx=sa, dx_by=sb, pred=pred)
        if key in ("diag", "axis_16"):
            Ns = E.cut_site(N, x0); ru = run(Ns, E57.sub(Bu, x0), "PS", 3.0, "xy"); rd = run(Ns, E57.sub(Bd, x0), "PS", 3.0, "xy")
            row["meas"] = {k: slope((np.array(ru[k]) - np.array(rd[k])) / 2, 0.1) for k in ("S", "L", "J")}; row["ratio"] = {k: row["meas"][k] / pred[k] for k in ("S", "J")}
        res[key] = row; print("60C", key, json.dumps(row), flush=True)
    return res


def part_D():
    n = 97; src, kap, N, Nd = E.source_terms(n); c = n // 2; x0 = (c + 16, c, c); Ns = E.cut_site(N, x0)
    Bu = E57.sub(E.shift(src["up"]["PS"], kap, "zero"), x0); Bd = E57.sub(E.shift(src["down"]["PS"], kap, "zero"), x0)
    ru = run(Ns, Bu, "PS", 3.0, "xy", dt=0.05); rd = run(Ns, Bd, "PS", 3.0, "xy", dt=0.05); om = slope((np.array(ru["S"]) - np.array(rd["S"])) / 2, 0.05)
    ref = json.loads((ROOT / "results/EX-GC59/EX-GC59_main.json").read_text())["runs"]["plane_16"]["Omega_meas"]
    res = dict(Omega_dt05=om, Omega_dt1=ref, rel=om / ref - 1); print("60D", res, flush=True); return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC60_{name}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    if job.startswith("A"):
        r = part_A(float(job[1:]))
    elif job == "B":
        r = part_B()
    elif job == "C":
        A3 = json.loads((OUT / "EX-GC60_A3.json").read_text())["rows"]; r = part_C(dict(rot=A3["xy_rot_PS"], shear=A3["xy_shear_PS"]))
    else:
        r = part_D()
    save(job, r, t0); print("done", job, "%.0fs" % (time.time() - t0), flush=True)
