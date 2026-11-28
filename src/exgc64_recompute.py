"""EX-GC64 recomputation -> results/EX-GC64/recompute.json: F1 by an explicit least-squares design matrix (full window and
t >= 2) from the saved momentum histories; precession by a separate lstsq on the saved angle histories. Shared: run records,
packet builders, operators, evolve routines, probe Hamiltonian. EXPLORATORY. Our construction; not LaViolette's Model G and not
Chain Theory."""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC64"; L = lambda k: json.loads((OUT / f"EX-GC64_{k}.json").read_text())
def F1(k, tmin=0.0):
    r = L(k)["runs"]; t = np.array(r["pair"]["t"]); y = np.array([a[0] - b[0] for a, b in zip(r["pair"]["Py"], r["one"]["Py"])]); m = t >= tmin - 1e-9
    return float(np.linalg.lstsq(np.stack([np.ones(m.sum()), t[m], t[m] ** 2], 1), y[m], rcond=None)[0][1])
res = {}
for tmin in (0.0, 2.0):
    g = lambda a, b, V: F1(a, tmin) / F1(b, tmin) / np.sqrt(1 - V * V); row = {}
    for tag in ("base", "shal", "deep", "narr"):
        row[tag] = dict(R_P=g(f"P_{tag}_2", f"P_{tag}_0", 0.2), R_E=g(f"E_{tag}_2", f"P_{tag}_0", 0.2))
    row["rev"] = dict(R_P=g("P_rev_2", "P_base_0", 0.2), R_E=g("E_rev_2", "P_base_0", 0.2)); row["slow"] = dict(R_P=g("P_slow_1", "P_base_0", 0.1), R_E=g("E_slow_1", "P_base_0", 0.1))
    row["RJE1"] = g("R_base_2", "P_base_0", 0.2); row["half"] = g("E_half_2", "E_half_0", 0.2); res[f"tmin_{tmin}"] = row
S = L("spin"); dt = 0.1; res["spin"] = {}
for k, v in S["runs"].items():
    om = {}
    for r in ("PS", "H"):
        a = (np.array(v[r]["phi_up"]) - np.array(v[r]["phi_down"])) / 2; t = dt * np.arange(len(a)); om[r] = float(np.linalg.lstsq(np.stack([np.ones_like(t), t], 1), a, rcond=None)[0][1])
    res["spin"][k] = om["H"] / om["PS"]
(OUT / "recompute.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
