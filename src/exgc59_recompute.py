"""EX-GC59 recomputation through separate paths -> results/EX-GC59/recompute.json.
(1) axis d = 12 and 16, PS +-S: probe rerun at dt = 0.05 with spin read from X.sdens sums (not vdot), slope by np.polyfit;
(2) the prediction from the LATTICE field instead of free-space GR: probe-density average of -1/2 (curl_op beta)_z of the solved shift;
(3) the light loop's coordinate loop integral -sum(bl)/(2A) from the solved shift, against the free-space predictor.
Shared: source construction, static solvers, probe builder. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G, exgc53 as T53, exgc57 as E57, exgc58 as X, exgc59 as E


def run(Ns, Bs, dt=0.05):
    psi = E.probe(); C = X.curl_op(Bs); H = lambda p: G.Hm(p, Ns, E.M, 0.0) + T53.Hb(p, Bs) + X.Hpl(p, C); ph = []; t = 0.0
    def rec(p):
        s = X.sdens(p); ph.append(np.arctan2(s[1].sum(), s[0].sum()))
    rec(psi)
    while t < E.T - 1e-9:
        k1 = -1j * H(psi); k2 = -1j * H(psi + dt / 2 * k1); k3 = -1j * H(psi + dt / 2 * k2); k4 = -1j * H(psi + dt * k3); psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; rec(psi)
    return np.array(ph), dt


if __name__ == "__main__":
    n = 97; src, kap, N, Nd = E.source_terms(n); Bu = E.shift(src["up"]["PS"], kap, "zero"); Bd = E.shift(src["down"]["PS"], kap, "zero"); res = {}
    main = json.loads((E.OUT / "EX-GC59_main.json").read_text())["runs"]; w = np.sum(np.abs(E.probe()) ** 2, -1); w /= w.sum(); Cg = X.curl_op(Bu)
    for key in ("axis_12", "axis_16"):
        x0 = E.pos(n, key); Ns = E.cut_site(N, x0); pu, dt = run(Ns, E57.sub(Bu, x0)); pd, _ = run(Ns, E57.sub(Bd, x0)); t = dt * np.arange(len(pu))
        om = float(np.polyfit(t, (pu - pd) / 2, 1)[0]); lat = float(np.sum(w * (-0.5) * E57.sub(Cg, x0)[2]))
        res[key] = dict(Omega_dt05=om, rel_to_main=om / main[key]["Omega_meas"] - 1, Omega_lattice_curl=lat, lattice_over_free=lat / main[key]["Omega_pred"], meas_over_lattice=om / lat)
        print(key, res[key], flush=True)
    x0 = E.pos(n, "axis_16"); st = T53.square_loop(x0, 8); _, bl = T53.loop_fields(st, N, Bu); A = 256.0
    res["light_loop"] = dict(Omega_loop_lattice=-float(np.sum(bl)) / (2 * A), free=main["axis_16"]["light"]["Omega_pred_loop"]); res["light_loop"]["ratio"] = res["light_loop"]["Omega_loop_lattice"] / res["light_loop"]["free"]
    print(res["light_loop"], flush=True); (E.OUT / "recompute.json").write_text(json.dumps(res, indent=1, default=float))
