"""EX-GC54 separate recomputation: R_dyn for every dynamic run from the raw P_y(t) lists, by an independent least-squares fit
(numpy lstsq on [1, t, t^2], and a second fit on t >= 2 only), compared with summary.json. Shares only the raw JSON records."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R54 = ROOT / "results" / "EX-GC54"
gam = 1 / np.sqrt(0.96)


def F1(rec, tmin=0.0):
    t = np.array(rec["runs"]["pair"]["t"]); P = np.array([p[0] for p in rec["runs"]["pair"]["Py"]]) - np.array([p[0] for p in rec["runs"]["one"]["Py"]])
    m = t >= tmin; A = np.stack([np.ones(m.sum()), t[m], t[m] ** 2], 1); c = np.linalg.lstsq(A, P[m], rcond=None)[0]
    return c[1]                     # dP/dt at t = 0 (as 53E)


S = json.loads((R54 / "summary.json").read_text()); out = {}
for r in ("D2q", "D2E"):
    a = json.loads((R54 / f"EX-GC54_B_{r}_0.0.json").read_text()); b = json.loads((R54 / f"EX-GC54_B_{r}_0.2.json").read_text())
    R = F1(b) / F1(a) * gam; Rlate = F1(b, 2.0) / F1(a, 2.0) * gam
    out[r] = dict(R=R, R_tge2=Rlate, summary=S["B"][r]["R_dyn"], diff=R - S["B"][r]["R_dyn"])
c = json.loads((R54 / "EX-GC54_C_P_0.2.json").read_text()); out["lam2"] = dict(F1=F1(c), summary=S["P5"]["F1"], diff_rel=F1(c) / S["P5"]["F1"] - 1)
(R54 / "recompute.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))
