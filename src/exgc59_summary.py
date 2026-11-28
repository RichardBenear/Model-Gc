"""EX-GC59 scoring -> results/EX-GC59/summary.json (registered P1-P7, C1-C4 of spec/EX-GC59_plan.md)."""
import json
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1] / "results" / "EX-GC59"; g = lambda n: json.loads((R / f"EX-GC59_{n}.json").read_text())
M = g("main"); r = M["runs"]; out = {}
out["P1"] = {k: r[k]["ratio"] for k in r}; out["P1"]["pass_"] = all(abs(v - 1) < 0.03 for v in out["P1"].values()) and r["axis_16"]["Omega_meas"] * r["plane_16"]["Omega_meas"] < 0
a = r["axis_16"]; ph = lambda k: np.array(a[k]["phi"]); z = ph("zero"); dp = (ph("PS_up") - z)[-1]
out["P2"] = dict(reversal=abs((ph("PS_up") - z)[-1] + (ph("PS_down") - z)[-1]) / abs(dp), unpol=abs((ph("unpol") - z)[-1]) / abs(dp))
out["P2"]["pass_"] = out["P2"]["reversal"] < 0.01 and out["P2"]["unpol"] < 0.01
out["P3"] = dict(P_over_PS=a["Omega_meas_P"] / a["Omega_meas"]); out["P3"]["pass_"] = abs(out["P3"]["P_over_PS"]) < 0.05
b129 = g("n129")["runs"]["axis_16"]; bd = g("dip")["runs"]["axis_16"]
out["P4"] = dict(n129_over_main=b129["Omega_meas"] / a["Omega_meas"], dip_over_main=bd["Omega_meas"] / a["Omega_meas"], n129_ratio=b129["ratio"], dip_ratio=bd["ratio"])
out["P4"]["pass_"] = abs(out["P4"]["n129_over_main"] - 1) < 0.02 and abs(out["P4"]["dip_over_main"] - 1) < 0.02 and abs(b129["ratio"] - 1) < 0.03 and abs(bd["ratio"] - 1) < 0.03
out["P5"] = dict(grad=abs((ph("grad") - z)[-1]) / abs(dp)); out["P5"]["pass_"] = out["P5"]["grad"] < 0.01
L = a["light"]; out["P6"] = dict(ratio=L["ratio"], sign_match=bool(np.sign(L["Omega_light"]) == np.sign(a["Omega_meas"])), note=L.get("note", ""))
out["P6"]["pass_"] = abs(L["ratio"] - 1) < 0.02 and out["P6"]["sign_match"]
As = g("async")["runs"]; om = [As[s]["Omega_meas"] for s in ("1", "2")]
out["P7"] = dict(async_over_ham=[o / a["Omega_meas"] for o in om], seeds_rel=abs(om[0] / om[1] - 1)); out["P7"]["pass_"] = all(abs(x - 1) < 1e-3 for x in out["P7"]["async_over_ham"]) and out["P7"]["seeds_rel"] < 1e-3
out["C1"] = dict(v=float(np.abs(ph("zero_Ndown") - z).max())); out["C1"]["pass_"] = out["C1"]["v"] < 1e-12
out["C2"] = dict(J=M["source"]["up"]["J_PS"], LS=M["source"]["up"]["L"] + M["source"]["up"]["Sz"]); out["C2"]["pass_"] = abs(out["C2"]["J"] / out["C2"]["LS"] - 1) < 1e-10
out["C3"] = dict(f=a["light_f"], pass_=abs(a["light_f"] / 1.045 - 1) < 0.01)
nh = max(a[k]["norm_dev"] for k in ("PS_up", "PS_down", "zero", "P_up")); na = max(As[s][x]["norm_dev"] for s in ("1", "2") for x in ("up", "down"))
out["C4"] = dict(ham=nh, events=na, pass_=nh < 1e-8 and na < 1e-12)
out["core_h"] = sum(g(n)["_meta"]["wall_s"] for n in ("main", "n129", "dip", "async")) / 3600
(R / "summary.json").write_text(json.dumps(out, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o)))
print({k: (v.get("pass_") if isinstance(v, dict) else v) for k, v in out.items()}); print(json.dumps(out, indent=0, default=float))
