"""EX-GC57 post-hoc diagnostics (not scored), run after the scored parts were seen:
(1) continuum linearised-GR prediction from the ring's actual momentum density (free space, no box):
    Omega_GR(x) = -(lam kappa / 8 pi) sum (x - r') x P(r') / |x - r'|^3, against the lattice point curl (P3's dipole -2 assumed a point source);
(2) the light clock's odd part in the ring sign, divided by the L = 64 instrument factor (C3), against the coordinate loop integral of its own path;
(3) the in-plane light loop at d = 24 with flat N, to see whether the part even in the ring sign comes from the N well;
(4) the orbiting probe against the coordinate loop integral of a staircase circle of its own radius (5), area by the shoelace formula.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc53 as T53
import exgc55_posthoc2 as SY
import exgc57 as E57
OUT = E57.OUT


def omega_free(P, kap, x, lam=4.0):
    n = P.shape[1]; LP = T53.link_pos(n); c = n // 2; xr = np.asarray(x, float) - c; tot = np.zeros(3)
    for a in range(3):
        m = np.abs(P[a]) > 1e-14 * np.abs(P).max()
        if not m.any():
            continue
        r = np.stack([LP[a][i][m] for i in range(3)], 1); d = xr - r; dd = np.linalg.norm(d, axis=1) ** 3
        v = np.zeros_like(r); v[:, a] = P[a][m]; tot += np.sum(np.cross(d, v) / dd[:, None], 0)
    return -(lam * kap / (8 * np.pi)) * tot


def staircase_area(steps):
    pts = [np.array(s[0], float) for s in steps]; x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


if __name__ == "__main__":
    n = 129; c = n // 2; q0, _ = T53.ring_sources(n, 0); kap = 4 * np.pi * 12 * 0.01 / float(q0.sum()); res = {}
    B57 = json.loads((OUT / "EX-GC57_B.json").read_text()); inst = B57["C3"]["uniform_over_2sum"]; F = {}
    for ell in (2, -2):
        q, P = T53.ring_sources(n, ell); J = T53.angmom(P); N, _ = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); F[ell] = (P, N, T53.shift_static(P, kap, J))
    P, N, B = F[2]
    for key, row in B57["positions"].items():
        kind, d = key.split("_"); d = int(d); x0 = (c, c, c + d) if kind == "axis" else (c + d, c, c); o = row["2"]; om = row["-2"]
        og = omega_free(P, kap, x0)[2]; r = dict(Omega_GR_free=og, curl_over_free=o["Omega_curl"][2] / og, orb_over_free=o["Omega_orb"] / og)
        dt_odd = (o["light"]["dt"] - om["light"]["dt"]) / 2; ro = o["light"]["readout"]
        r["light_odd_over_readout"] = dt_odd / ro; r["light_odd_over_readout_inst"] = dt_odd / ro / inst
        st = T53.staircase_loop(x0, 5.0, (0, 0, 1)); Ns, bl = T53.loop_fields(st, N, B); A = staircase_area(st)
        r["Omega_loop_r5"] = -float(np.sum(bl)) / (2 * A); r["orb_over_loop_r5"] = o["Omega_orb"] / r["Omega_loop_r5"]; r["area_r5"] = A
        res[key] = r; print(key, json.dumps(r, default=float), flush=True)
    res["ratio_axis32_plane32_free"] = res["axis_32"]["Omega_GR_free"] / res["plane_32"]["Omega_GR_free"]
    x0 = (c + 24, c, c); one = np.ones((n, n, n)); fl = {}
    for ell in (2, -2):
        st = T53.square_loop(x0, 8); Ns, bl = T53.loop_fields(st, one, F[ell][2]); fl[str(ell)] = SY.wave(Ns, bl, +1, -1, width=4.0) - SY.wave(Ns, bl, -1, -1, width=4.0)
    fl["even_over_odd"] = (fl["2"] + fl["-2"]) / (fl["2"] - fl["-2"]); res["plane24_flatN"] = fl
    print("ratio free", res["ratio_axis32_plane32_free"], "flatN", fl, flush=True)
    (OUT / "posthoc.json").write_text(json.dumps(res, indent=1, default=float))
