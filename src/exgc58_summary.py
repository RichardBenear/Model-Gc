"""EX-GC58 scoring -> results/EX-GC58/summary.json (registered P1-P6, C1-C4 of spec/EX-GC58_plan.md)."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results" / "EX-GC58"; g = lambda n: json.loads((R / f"EX-GC58_{n}.json").read_text())
gam = lambda V: 1 / np.sqrt(1 - V * V); out = {}; names = []
A = g("A"); names.append("A")
nd = max(r[m]["norm_dev"] for r in A["rows"] for m in ("P", "PS"))
out["P1"] = dict(norm_dev=nd, slopes=A["slopes"], pass_=nd < 1e-13 and all(abs(v - 2) < 0.15 for s in A["slopes"].values() for v in s.values()))
Bs = [g(f"B{i}") for i in range(10)]; names += [f"B{i}" for i in range(10)]
out["P2"] = dict(cases=[{k: b[k] for k in ("m", "k0", "r", "chi", "ratio", "rel_ham")} for b in Bs])
out["P2"]["max_dev_pred"] = max(abs(b["ratio"] - 1) for b in Bs); out["P2"]["max_rel_ham"] = max(abs(b["rel_ham"]) for b in Bs)
out["P2"]["pass_"] = out["P2"]["max_dev_pred"] < 2e-3 and out["P2"]["max_rel_ham"] < 1e-4
C = g("C"); names.append("C"); pos = C["positions"]; ax = ["axis_24", "axis_32"]
ps = [pos[k]["2"][f"spin_{s}_PS_dEdJ"] for k in ax for s in ("3", "5")]; pp = [pos[k]["2"][f"spin_{s}_P_dEdJ"] for k in ax for s in ("3", "5")]
orb = [abs(pos[k]["2"]["orb_PS"] / pos[k]["2"]["orb_P"] - 1) for k in pos]
out["P3"] = dict(PS_dEdJ=ps, P_dEdJ=pp, orb_change_max=max(orb), grad_rel=C["grad"]["rel"])
out["P3"]["pass_"] = all(abs(x - 1) < 0.03 for x in ps) and all(x < 0.15 for x in pp) and max(orb) < 0.01 and C["grad"]["rel"] < 1e-12
D = g("D"); names.append("D")
LS = D["L"][2] + D["S"][2]; dl = [D[r]["loops"][d]["ratio"] for r in ("P", "PS") for d in ("20", "28")]
out["P4"] = dict(J_PS_over_LS=D["PS"]["J"][2] / LS, J_P_over_L=D["P"]["J"][2] / D["L"][2], loop_over_dipole=dl, J_PS=D["PS"]["J"][2])
out["P4"]["pass_"] = abs(out["P4"]["J_PS_over_LS"] - 1) < 0.01 and abs(out["P4"]["J_P_over_L"] - 1) < 0.01 and all(abs(x - 1) < 0.05 for x in dl)
E1, E2 = g("E_0.1"), g("E_0.05"); names += ["E_0.1", "E_0.05"]
out["P5"] = dict(drift_rel=E1["drift_rel"], drift_rel_half=E2["drift_rel"], fall=E1["drift"] / E2["drift"])
out["P5"]["pass_"] = E1["drift_rel"] < 1e-4 and out["P5"]["fall"] >= 10
F0, F2 = g("F_PS_0.0"), g("F_PS_0.2"); names += ["F_PS_0.0", "F_PS_0.2"]
C55 = {k: json.loads((ROOT / f"results/EX-GC55/EX-GC55_C_P_{k}.json").read_text())["F1"] for k in ("0.0", "0.2")}
Rps = F2["F1"] / F0["F1"] * gam(0.2); out["P6"] = dict(R_PS=Rps, R_P55=C55["0.2"] / C55["0.0"] * gam(0.2), dF1_V0=F0["F1"] / C55["0.0"] - 1, dF1_V02=F2["F1"] / C55["0.2"] - 1)
out["P6"]["pass_"] = abs(Rps - 1) < 0.015 and abs(Rps - 0.9986) < 0.005
out["C1"] = dict(v=A["C1"], pass_=A["C1"] < 1e-15); out["C2"] = dict(A["C2"], pass_=A["C2"]["herm"] < 1e-12 and A["C2"]["rel"] < 1e-12)
out["C3"] = dict(d24=[pos["axis_24"]["2"]["spin_3_P_57def"], pos["axis_24"]["2"]["spin_5_P_57def"]], ref57=[0.10942915009187247, 0.04775967089796942])
out["C3"]["pass_"] = all(abs(a / b - 1) < 1e-8 for a, b in zip(out["C3"]["d24"], out["C3"]["ref57"]))
try:
    FP = g("F_P_0.2"); names.append("F_P_0.2"); out["C4"] = dict(F1=FP["F1"], F1_55=C55["0.2"], rel=abs(FP["F1"] / C55["0.2"] - 1)); out["C4"]["pass_"] = out["C4"]["rel"] < 1e-10
except FileNotFoundError:
    out["C4"] = "skipped"
out["core_h"] = sum(g(n)["_meta"]["wall_s"] for n in names) / 3600
(R / "summary.json").write_text(json.dumps(out, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o))); print({k: (v.get("pass_") if isinstance(v, dict) else v) for k, v in out.items()})
