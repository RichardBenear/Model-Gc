"""EX-GC55 separate recomputation from the raw run records: R_dyn by an independent lstsq fit of P1(t) (and t >= 2 only),
and the energy drift from the stored E(t) components (Em + Ei + Eb re-summed). Compared with summary.json."""
import json
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1] / "results" / "EX-GC55"; g = lambda n: json.loads((R / f"EX-GC55_{n}.json").read_text())
gam = 1 / np.sqrt(0.96); S = json.loads((R / "summary.json").read_text()); out = {}
def F1(rec, tmin=0.0):
    t = np.array(rec["runs"]["pair"]["t"]); P = np.array([p[0] for p in rec["runs"]["pair"]["Py"]]) - np.array([p[0] for p in rec["runs"]["one"]["Py"]])
    m = t >= tmin; A = np.stack([np.ones(m.sum()), t[m], t[m] ** 2], 1); return np.linalg.lstsq(A, P[m], rcond=None)[0][1]
for r in ("P", "RJE"):
    a, b = g(f"C_{r}_0.0"), g(f"C_{r}_0.2")
    out[r] = dict(R=F1(b) / F1(a) * gam, R_tge2=F1(b, 2) / F1(a, 2) * gam, summary=S["P4"][f"R_{r}"])
out["diff"] = out["RJE"]["R"] - out["P"]["R"]; out["diff_tge2"] = out["RJE"]["R_tge2"] - out["P"]["R_tge2"]
for r in ("P", "MIX", "RJE", "P_half"):
    h = g(f"B_{r}")["hist"]; E = np.array([x["Em"] + x["Ei"] + x["Eb"] for x in h]); out[f"drift_{r}"] = float(np.max(np.abs(E - E[0])) / abs(h[0]["Ei"]))
(R / "recompute.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))
