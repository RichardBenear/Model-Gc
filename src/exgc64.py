"""EX-GC64: transfer, reciprocity and hybrid spin for 63's per-local-tick count (D2-E1). 64A/B pair runs (55 box: n = 49,
sep 14, T = 8) for D2-E1 and P at several depth / width / boost / dt settings; 64C reciprocal flat energy current (RJE1);
64D 59 probes with the hybrid source s(D2-E1) + 1/4 C^T s; C1-C4. Plan: spec/EX-GC64_plan.md. EXPLORATORY. Our construction;
not LaViolette's Model G and not Chain Theory."""
import json, sys, time, platform, math
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G, exgc52 as S, exgc53 as T53, exgc54 as X54, exgc55 as X55, exgc63 as Q63
import exgc59 as E59, exgc57 as E57, exgc58 as X

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC64"; Nf = G.Nf; LAM = 4.0
_c55 = X55.coupling


def je1_terms(psi, m, B):
    return X55.je_terms(psi, np.ones(psi.shape[:3]), m, B)


def coupling(psi, N, m, B, rule):
    if rule == "RJE1":
        Hb, s = je1_terms(psi, m, B); return Hb, s, s
    return _c55(psi, N, m, B, rule)


X55.coupling = coupling        # 55's evolve looks the name up at call time


def setup(V, kscale=1.0, sigma=None, n=49, sep=14):
    c = n // 2; m = S.M; sig = S.SIG if sigma is None else sigma; k0 = np.sign(V) * S.k_for_v(abs(V), m, 0.0)
    p1, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), sig, (k0, 0, 0)); p2, _, _ = S.packet_k(n, m, 0.0, (c, c + sep // 2, c), sig, (k0, 0, 0))
    p0, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), sig, (0, 0, 0)); kap = kscale * 4 * np.pi * sep * 0.01 / float(G.charge(p0, np.ones((n, n, n)), m, 0.0).sum())
    return p1, p2, kap, m, c


def strength(V, kscale=1.0, sigma=None, n=49):
    p1, p2, kap, m, c = setup(V, kscale, sigma, n); ps = [p1, p2]; f, g, info = S.initial_field(ps, [m, m], 0.0, kap, n, c, V); N = Nf(f)
    P = sum(T53.link_P(p) for p in ps)[0].sum()
    return dict(D2E1=float(4 * Q63.source(ps, N, [m, m], "D2E1")[0].sum() / P), D2E=float(4 * X54.source(ps, N, [m, m], "D2E")[0].sum() / P),
                RJE1_src=float(4 * sum(je1_terms(p, m, np.zeros((3, n, n, n)))[1] for p in ps)[0].sum() / P), Nmin=float(N.min()))


def run_pair(V, rule, kscale=1.0, sigma=None, dt=0.1, n=49, T=8.0):
    p1, p2, kap, m, c = setup(V, kscale, sigma, n); out = dict(V=V, rule=rule, kscale=kscale, sigma=sigma, dt=dt, n=n, T=T, kappa=kap, runs={})
    for lab, ps in (("pair", [p1, p2]), ("one", [p1])):
        t0 = time.time(); f, g, info = S.initial_field(ps, [m] * len(ps), 0.0, kap, n, c, V); N = Nf(f); z = np.zeros((3, n, n, n))
        if rule == "RJE1":
            j0 = sum(coupling(p, N, m, z, rule)[1] for p in ps); B, gB = X54.initial_shift(j0, kap, LAM, n, V)
            _, hist = X55.evolve(ps, [m] * len(ps), f, g, B, gB, kap, LAM, T, dt, rule); hist = [(h["t"], h["Py"]) for h in hist]
        else:
            j0 = X54.source(ps, N, [m] * len(ps), rule); B, gB = X54.initial_shift(j0, kap, LAM, n, V)
            _, hist = X54.evolve(ps, [m] * len(ps), f, g, B, gB, kap, LAM, T, dt, rule)
        out["runs"][lab] = dict(t=[h[0] for h in hist], Py=[h[1] for h in hist], wall=time.time() - t0, j0_tot=j0.sum(axis=(1, 2, 3)).tolist())
        print("64 %s V %.2f k %.2f sig %s dt %.3f %s %.0fs" % (rule, V, kscale, sigma, dt, lab, time.time() - t0), flush=True)
    t = np.array(out["runs"]["pair"]["t"]); P1 = np.array([p[0] for p in out["runs"]["pair"]["Py"]]) - np.array([p[0] for p in out["runs"]["one"]["Py"]])
    out["F1"] = float(np.polyfit(t, P1, 2)[1]); return out


def part_spin(n=97, keys=("axis_12", "axis_16", "plane_16")):
    src, kap, N, Nd = E59.source_terms(n); res = dict(kappa=kap, runs={}, mz={})
    for s, chi in (("up", E59.UP), ("down", E59.DOWN)):
        p = E59.source(n, chi); src[s]["H"] = src[s]["PS"] - src[s]["P"] + Q63.d2e1(p, E59.M)
        res["mz"][s] = {r: float(Q63.moment(src[s][r], n)[2]) for r in ("PS", "H", "P")}
    B = {(r, s): E59.shift(src[s][r], kap, "zero") for r in ("PS", "H") for s in ("up", "down")}
    for key in keys:
        x0 = E59.pos(n, key); Ns = E59.cut_site(N, x0); row = {}
        for r in ("PS", "H"):
            ph = {s: np.array(E59.run_ham(Ns, E57.sub(B[(r, s)], x0), "PS")["phi"]) for s in ("up", "down")}
            row[r] = dict(Omega=E59.slope((ph["up"] - ph["down"]) / 2), up=E59.slope(ph["up"]), down=E59.slope(ph["down"]), phi_up=ph["up"].tolist(), phi_down=ph["down"].tolist())
        row["ratio"] = row["H"]["Omega"] / row["PS"]["Omega"]; res["runs"][key] = row; print("64D", key, "ratio %.5f" % row["ratio"], flush=True)
    return res


def part_checks(n=49):
    res = {}; m = S.M; c = n // 2; S_ = S
    p1, p2, kap, _, _ = setup(0.2, n=n); k0 = S.k_for_v(0.2, m, 0.0)
    q1, _, _ = S.packet_k(n, m, 0.0, (c, c - 7, c), S.SIG, (-k0, 0, 0)); s1 = Q63.d2e1(p1, m); sq = Q63.d2e1(q1, m)
    a = [math.fsum(s1[0].ravel()), math.fsum(sq[0].ravel())]; res["C1_fsum_flip"] = abs(a[0] + a[1]) / abs(a[0])
    Pm = np.array(G.AL[1]) @ np.array(G.AL[2]) @ np.array(G.BE)       # parity: anticommutes with alpha_x, commutes with alpha_y, alpha_z, beta
    pr = np.einsum("ab,...b->...a", Pm, p1[::-1].copy())             # x -> (n-1) - x
    sr = Q63.d2e1(pr, m); sx_ref = -np.roll(s1[0][::-1], -1, axis=0)  # link (x, x+1) maps to link (n-2-x, n-1-x)
    res["C1_parity_local_x"] = float(np.abs(sr[0][:-1] - sx_ref[:-1]).max() / np.abs(s1[0]).max())
    res["C1_parity_local_yz"] = float(max(np.abs(sr[a][:] - s1[a][::-1]).max() for a in (1, 2)) / np.abs(s1[0]).max())
    res["C1_overlap_q1_vs_parity"] = float(abs(np.vdot(q1, pr)) / (np.linalg.norm(q1) * np.linalg.norm(pr)))
    rng = np.random.default_rng(64); N1 = np.ones((13,) * 3); B = 1e-3 * rng.standard_normal((3, 13, 13, 13)); ps = rng.standard_normal((13, 13, 13, 4)) + 1j * rng.standard_normal((13, 13, 13, 4))
    h1, s_1 = je1_terms(ps, m, B); h0, s_0 = X55.je_terms(ps, N1, m, B)
    res["C3_RJE1_vs_RJE"] = float(max(np.abs(h1 - h0).max() / np.abs(h0).max(), np.abs(s_1 - s_0).max() / np.abs(s_0).max()))
    res["C3_D2E1_vs_MIX"] = float(np.abs(Q63.d2e1(ps, m) - _c55(ps, N1, m, B, "MIX")[1]).max() / np.abs(Q63.d2e1(ps, m)).max())
    res["strength"] = {lab: strength(0.2, kscale=k, sigma=sg) for lab, k, sg in (("base", 1, None), ("shallow", 0.5, None), ("deep", 2, None), ("narrow", 1, 4.0))}
    print(json.dumps(res, indent=1), flush=True); return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC64_{name}.json").write_text(json.dumps(res, indent=1, default=float))


JOBS = {  # name: (V, rule, kscale, sigma, dt)
    "P_base_0": (0.0, "P", 1, None, 0.1), "P_base_2": (0.2, "P", 1, None, 0.1), "E_base_2": (0.2, "D2E1", 1, None, 0.1), "E_base_0": (0.0, "D2E1", 1, None, 0.1),
    "P_shal_0": (0.0, "P", 0.5, None, 0.1), "P_shal_2": (0.2, "P", 0.5, None, 0.1), "E_shal_2": (0.2, "D2E1", 0.5, None, 0.1),
    "P_deep_0": (0.0, "P", 2, None, 0.1), "P_deep_2": (0.2, "P", 2, None, 0.1), "E_deep_2": (0.2, "D2E1", 2, None, 0.1),
    "P_narr_0": (0.0, "P", 1, 4.0, 0.1), "P_narr_2": (0.2, "P", 1, 4.0, 0.1), "E_narr_2": (0.2, "D2E1", 1, 4.0, 0.1),
    "P_slow_1": (0.1, "P", 1, None, 0.1), "E_slow_1": (0.1, "D2E1", 1, None, 0.1),
    "P_rev_2": (-0.2, "P", 1, None, 0.1), "E_rev_2": (-0.2, "D2E1", 1, None, 0.1),
    "R_base_2": (0.2, "RJE1", 1, None, 0.1),
    "E_half_0": (0.0, "D2E1", 1, None, 0.05), "E_half_2": (0.2, "D2E1", 1, None, 0.05)}

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    if job == "checks":
        save(job, part_checks(), t0)
    elif job == "spin":
        save(job, part_spin(), t0)
    else:
        V, rule, k, sg, dt = JOBS[job]; save(job, run_pair(V, rule, k, sg, dt), t0)
