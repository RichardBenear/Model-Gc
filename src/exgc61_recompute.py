"""EX-GC61 recomputation -> results/EX-GC61/recompute.json.
(1) P5/P6: axis-16 and diagonal probes rerun for one source spin per call (job up | down) with field and probe at dt = 0.05, spin read
    from X.sdens sums, J from sdens + link_P, odd rate by np.polyfit; (2) arrival times from a cubic-spline crossing of the saved site
    curls (job arr). Shared: source builder, field class, probe Hamiltonian, frozen coefficients. EXPLORATORY. Our construction; not
    LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc61 as Q, exgc59 as E, exgc57 as E57, exgc58 as X, exgc53 as T53, exgc60 as P60


class Probe2(Q.Probe):
    def add(self):
        s = X.sdens(self.psi).sum(axis=(1, 2, 3)); L = np.array(T53.angmom(T53.link_P(self.psi)), float); self.rec.append(np.concatenate([s, L, L + 0.5 * s]).tolist())


def probes(spin, dt=0.05):
    src, kap, N, Nd = E.source_terms(97); F = Q.Field(src[spin]["PS"], Q.LAM * kap); c = 48; xs = {"axis_16": (c, c, c + 16), "diag": (c + 12, c, c + 12)}; live = {}
    while F.t < 50.0 - 1e-9:
        if abs(F.t - 40.0) < 1e-6:
            live = {k: Probe2(E.cut_site(N, x), P60.probe(3.0, "x")) for k, x in xs.items()}
        olds = {k: E57.sub(F.B, xs[k]) for k in live}; F.step(dt)
        for k, pr in live.items():
            pr.step(olds[k], E57.sub(F.B, xs[k]), dt)
    return {k: pr.rec for k, pr in live.items()}


if __name__ == "__main__":
    job = sys.argv[1]; out = Q.OUT
    if job in ("up", "down"):
        (Q.PART / f"rc_{job}.json").write_text(json.dumps(probes(job)))
    else:
        from scipy.interpolate import CubicSpline
        B = json.loads((out / "EX-GC61_B97.json").read_text()); res = {"t_half_spline": {}}
        for k, om in B["om"].items():
            om = np.array(om); t = Q.DT * np.arange(len(om)); ref = om[500]; cs = CubicSpline(t, om / ref - 0.5); r = [x for x in cs.roots() if 0 < x < 50]
            res["t_half_spline"][k] = float(r[0])
        U = json.loads((Q.PART / "rc_up.json").read_text()); D = json.loads((Q.PART / "rc_down.json").read_text()); S = json.loads((out / "summary.json").read_text())
        for k in ("axis_16", "diag"):
            u = np.array(U[k]); d = np.array(D[k]); t = 0.05 * np.arange(len(u)); row = {}
            for nm, off in (("S", 0), ("J", 6)):
                a = (np.arctan2(u[:, off + 1], u[:, off]) - np.arctan2(d[:, off + 1], d[:, off])) / 2; row[nm] = float(np.polyfit(t, a, 1)[0])
            main = S["P5"]["rate"] if k == "axis_16" else S["P6"]["rate"]; row["rel_to_main"] = {nm: row[nm] / main[nm] - 1 for nm in ("S", "J")}; res[k] = row
        res["t_half_rel_diff_max"] = max(abs(res["t_half_spline"][k] - B["t_half"][k]) for k in B["t_half"])
        print(res); (out / "recompute.json").write_text(json.dumps(res, indent=1, default=float))
