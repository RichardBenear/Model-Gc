"""EX-GC56 separate recomputation: R_dyn of the 56C runs from raw P_y(t) by independent lstsq (also t >= 2), and the 56A
energy-current band prediction from the closed form 1/2 d^2(E^2)/dk^2 = cos 2k + r^2 sin^2 k + (m + r(1 - cos k)) r cos k (k-space
weights rebuilt), compared with the stored finite-difference value. Shares only exgc52.packet_k and the run records."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc52 as S
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results" / "EX-GC56"; g = lambda n: json.loads((R / f"EX-GC56_{n}.json").read_text())
gam = lambda V: 1 / np.sqrt(1 - V * V); out = {}
def F1(rec, tmin=0.0):
    t = np.array(rec["runs"]["pair"]["t"]); P = np.array([p[0] for p in rec["runs"]["pair"]["Py"]]) - np.array([p[0] for p in rec["runs"]["one"]["Py"]])
    m = t >= tmin; A = np.stack([np.ones(m.sum()), t[m], t[m] ** 2], 1); return np.linalg.lstsq(A, P[m], rcond=None)[0][1]
c55 = lambda k: json.loads((ROOT / f"results/EX-GC55/EX-GC55_C_P_{k}.json").read_text())
for tm in (0.0, 2.0):
    out[f"tmin{tm}"] = dict(R_kappa_half=F1(g("Ck_0.2"), tm) / F1(g("Ck_0.0"), tm) * gam(0.2), R_V015=F1(g("Cv_0.15"), tm) / F1(c55("0.0"), tm) * gam(0.15),
                            R_sigma4=F1(g("Cs_0.2"), tm) / F1(g("Cs_0.0"), tm) * gam(0.2), R55=F1(c55("0.2"), tm) / F1(c55("0.0"), tm) * gam(0.2))
rows = []
for x in [y for y in g("A")["matter"] if y["rule"] == "JE" and y["chi"][0] == 1]:
    _, K, w = S.packet_k(48, x["m"], x["r"], (24, 24, 24), S.SIG, (x["k0"], 0, 0)); k = K[..., 0]; r = x["r"]; m = x["m"]
    cf = -float(np.sum(w * (np.cos(2 * k) + r * r * np.sin(k) ** 2 + (m + r * np.sum(1 - np.cos(K), -1)) * r * np.cos(k))))
    rows.append(dict(m=m, r=r, k0=x["k0"], closed=cf, stored=x["pred_JE"], rel=cf / x["pred_JE"] - 1))
out["JE_closed_form"] = rows; out["JE_closed_max_rel"] = max(abs(r_["rel"]) for r_ in rows)
(R / "recompute.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps({k: v for k, v in out.items() if k != "JE_closed_form"}, indent=1, default=float))
