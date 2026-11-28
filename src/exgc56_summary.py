"""EX-GC56 scoring -> results/EX-GC56/summary.json (registered P1-P5, C1-C3 of spec/EX-GC56_plan.md)."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results" / "EX-GC56"; g = lambda n: json.loads((R / f"EX-GC56_{n}.json").read_text())
gam = lambda V: 1 / np.sqrt(1 - V * V); out = {}
A = g("A"); M = A["matter"]
P = [x for x in M if x["rule"] == "P"]; J = [x for x in M if x["rule"] == "JE"]
rp = np.array([x["ratio_P"] for x in P])
out["P1"] = dict(ratios=rp.tolist(), max_dev=float(np.max(np.abs(rp - 1))), spread=float(rp.max() - rp.min()))
out["P1"]["pass_"] = out["P1"]["max_dev"] < 2e-3 and out["P1"]["spread"] < 2e-3
L = A["light"]; out["P2"] = dict(ratios=[x["ratio"] for x in L], pass_=all(abs(x["ratio"] - 1) < 5e-3 for x in L))
rj = np.array([x["ratio_JE"] for x in J]); key = [x for x in J if x["r"] == 1.0 and x["m"] == 0.6 and x["k0"] == 0.0][0]
keyP = [x for x in P if x["r"] == 1.0 and x["m"] == 0.6 and x["k0"] == 0.0][0]
out["P3"] = dict(ratios_JE_band=rj.tolist(), max_dev_band=float(np.max(np.abs(rj - 1))), key_JE_over_P=key["drift"] / keyP["drift"],
                 JE_over_P_all=[x["drift"] / y["drift"] for x, y in zip(J, P)])
out["P3"]["pass_"] = out["P3"]["max_dev_band"] < 0.01 and abs(out["P3"]["key_JE_over_P"] - 1) > 0.30
B = g("B"); out["P4"] = dict(P_identity=B["P_identity"], JE_residual=B["JE_residual"], pass_=B["P_identity"] < 1e-14 and B["JE_residual"] > 1e-2)
F = {k: g(k)["F1"] for k in ("Ck_0.0", "Ck_0.2", "Cv_0.15", "Cs_0.0", "Cs_0.2", "C3_0.2")}
C55 = {k: json.loads((ROOT / f"results/EX-GC55/EX-GC55_C_P_{k}.json").read_text())["F1"] for k in ("0.0", "0.2")}
Rk = F["Ck_0.2"] / F["Ck_0.0"] * gam(0.2); Rv = F["Cv_0.15"] / C55["0.0"] * gam(0.15); Rs = F["Cs_0.2"] / F["Cs_0.0"] * gam(0.2); R55 = C55["0.2"] / C55["0.0"] * gam(0.2)
out["P5"] = dict(R_kappa_half=Rk, R_V015=Rv, R_sigma4=Rs, R55=R55, dk=Rk - R55, ds=Rs - R55)
out["P5"]["pass_"] = all(abs(x - 1) < 0.01 for x in (Rk, Rv, Rs)) and abs(Rk - R55) < 0.005 and abs(Rs - R55) < 0.01
out["C1"] = A["C1"]; out["C2"] = A["C2"]; out["C3"] = dict(F1=F["C3_0.2"], F1_55=C55["0.2"], rel=abs(F["C3_0.2"] / C55["0.2"] - 1))
out["C3"]["pass_"] = out["C3"]["rel"] < 1e-10
out["core_h"] = sum(g(n)["_meta"]["wall_s"] for n in ["A", "B", "Ck_0.0", "Ck_0.2", "Cv_0.15", "Cs_0.0", "Cs_0.2", "C3_0.2"]) / 3600
(R / "summary.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))
