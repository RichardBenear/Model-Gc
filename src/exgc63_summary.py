"""EX-GC63 scoring -> results/EX-GC63/summary.json, gates exactly as spec/EX-GC63_plan.md. EXPLORATORY. Our construction; not
LaViolette's Model G and not Chain Theory."""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC63"; L = lambda k: json.loads((OUT / f"EX-GC63_{k}.json").read_text())
def F1(V):
    p = L(f"A_{V}_pair")["runs"]["pair"]; o = L(f"A_{V}_one")["runs"]["one"]; t = np.array(p["t"])
    return float(np.polyfit(t, np.array([x[0] for x in p["Py"]]) - np.array([x[0] for x in o["Py"]]), 2)[1])
C = L("checks"); B = L("B"); f0, f2 = F1("0.0"), F1("0.2"); R = f2 / f0 / np.sqrt(1 - 0.04); res = {}
res["P1"] = dict(v=C["P1_strength_pair"], pass_=bool(abs(C["P1_strength_pair"] - 3.85) <= 0.02))
res["P2"] = dict(R=R, F1_0=f0, F1_02=f2, line_pred=1.1152 - 0.02943 * C["P1_strength_pair"], pass_=bool(abs(R - 1.002) <= 0.008))
res["P3"] = dict(v=abs(R - 1), pass_=bool(abs(R - 1) < 0.02))
res["P4"] = dict(rel=f0 / 4.699e-4 - 1, pass_=bool(abs(f0 / 4.699e-4 - 1) < 0.01))
res["P5"] = dict(up=B["up"]["ratio"], down=B["down"]["ratio"], pass_=bool(B["up"]["ratio"] < 0.05 and B["down"]["ratio"] < 0.05))
res["C1"] = dict(v=C["C1_zero"], pass_=bool(C["C1_zero"] == 0)); res["C2"] = dict(v=C["C2_flip"], pass_=bool(C["C2_flip"] < 1e-12))
res["C3"] = dict(v=C["C3_flatN"], pass_=bool(C["C3_flatN"] < 1e-14))
res["core_h"] = sum(L(k)["_meta"]["wall_s"] for k in ("A_0.0_pair", "A_0.2_pair", "A_0.0_one", "A_0.2_one", "B", "checks")) / 3600
(OUT / "summary.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
