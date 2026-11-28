"""EX-GC58 recomputation through separate paths -> results/EX-GC58/recompute.json.
P1: the vectorised class sweep (not the sequential link loop) against the exact exponential, slopes over s.
P3: spin energies as operator expectations <psi|(H_beta + H_pl) psi> with the curl built on the probe box (not cut from the global curl).
P4: S from vdot(psi, Sigma_z psi); the dipole J recomputed as angmom(P) + angmom(1/4 C^T s) separately.
P5/P6: drift from the stored energy histories; F1 by an independent lstsq on the raw momenta (all t and t >= 2).
Shared dependencies: 53's ring/probe/shift code, 55's evolve. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc57 as E57
import exgc58 as X
R = X.OUT; ROOT = X.ROOT; g = lambda n: json.loads((R / f"EX-GC58_{n}.json").read_text())


def f1(run, tmin=0.0):
    t = np.array(run["runs"]["pair"]["t"]); P1 = np.array([p[0] for p in run["runs"]["pair"]["Py"]]) - np.array([p[0] for p in run["runs"]["one"]["Py"]])
    m = t >= tmin; A = np.stack([np.ones(m.sum()), t[m], t[m] ** 2], 1); return float(np.linalg.lstsq(A, P1[m], rcond=None)[0][1])


if __name__ == "__main__":
    res = {}; rng = np.random.default_rng(58); n = 13; dt = 0.1
    psi0 = rng.standard_normal((n, n, n, 4)) + 1j * rng.standard_normal((n, n, n, 4)); psi0 /= np.linalg.norm(psi0); B1 = rng.standard_normal((3, n, n, n)) * 1e-3
    errs = []
    for s in (1.0, 2.0, 4.0):
        B = s * B1; ex = X.expm_apply(lambda p: T53.Hb(p, B), psi0, dt); errs.append(float(np.linalg.norm(X.class_sweep(psi0.copy(), B, dt) - ex)))
    res["P1_class_sweep_slope"] = float(np.polyfit(np.log([1, 2, 4]), np.log(errs), 1)[0]); print("P1", res["P1_class_sweep_slope"], flush=True)
    nn = 129; c = nn // 2; q0, _ = T53.ring_sources(nn, 0); kap = 4 * np.pi * 12 * 0.01 / float(q0.sum()); q, P = T53.ring_sources(nn, 2); B = T53.shift_static(P, kap, T53.angmom(P))
    Sz = lambda p: float(np.real(np.vdot(p, p @ X.SIG[2].T))); res["P3"] = {}
    for d in (24, 32):
        x0 = (c, c, c + d); Bs = E57.sub(B, x0); Cs = X.curl_op(Bs); oc = -0.5 * E57.curl_site(B, x0)[2]
        for sig in (3.0, 5.0):
            up, dn = E57.spin_packet(sig, T53.UP), E57.spin_packet(sig, T53.DOWN); H = lambda p: float(np.real(np.vdot(p, T53.Hb(p, Bs) + X.Hpl(p, Cs))))
            J = lambda p: float(T53.angmom(T53.link_P(p))[2]) + 0.5 * Sz(p)
            res["P3"][f"{d}_{sig:g}"] = (H(up) - H(dn)) / (J(up) - J(dn)) / oc; print("P3", d, sig, res["P3"][f"{d}_{sig:g}"], flush=True)
    nd = 97; cd = nd // 2; psi, _, _ = S.packet_k(nd, S.M, 0.0, (cd, cd, cd), 3.0, (0, 0, 0), chi=T53.UP); Lc = T53.angmom(T53.link_P(psi))[2]; Jsp = T53.angmom(0.25 * X.curl_T(X.sdens(psi)))[2]
    res["P4"] = dict(L=Lc, J_spin_part=Jsp, half_Sz=0.5 * Sz(psi), ratio=Jsp / (0.5 * Sz(psi))); print("P4", res["P4"], flush=True)
    E = [np.array([h["E"] for h in g(k)["hist"]]) for k in ("E_0.1", "E_0.05")]; Ei = abs(g("E_0.1")["hist"][0]["Ei"])
    res["P5"] = dict(drift_rel=float(np.max(np.abs(E[0] - E[0][0])) / Ei), fall=float(np.max(np.abs(E[0] - E[0][0])) / np.max(np.abs(E[1] - E[1][0]))))
    gam = 1 / np.sqrt(1 - 0.04); F0, F2 = g("F_PS_0.0"), g("F_PS_0.2")
    res["P6"] = dict(R=f1(F2) / f1(F0) * gam, R_t2=f1(F2, 2.0) / f1(F0, 2.0) * gam); print("P5", res["P5"], "P6", res["P6"], flush=True)
    (R / "recompute.json").write_text(json.dumps(res, indent=1, default=float))
