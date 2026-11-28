"""EX-GC64 scoring -> results/EX-GC64/summary.json, gates exactly as spec/EX-GC64_plan.md. EXPLORATORY. Our construction; not
LaViolette's Model G and not Chain Theory."""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC64"; L = lambda k: json.loads((OUT / f"EX-GC64_{k}.json").read_text())
C = L("checks"); st = C["strength"]; K = 0.02943; res = {}; F = lambda k: L(k)["F1"]
R = lambda k2, k0, V: F(k2) / F(k0) / np.sqrt(1 - V * V)
res["P1"] = dict(s={k: v["D2E1"] for k, v in st.items()}, D2E_54={k: v["D2E"] for k, v in st.items()},
                 pass_=bool(all(abs(st[k]["D2E1"] - 3.85) <= 0.02 for k in ("base", "shallow", "deep"))))
rows = {}
for lab, tag, sk, V in (("base", "base", "base", 0.2), ("shallow", "shal", "shallow", 0.2), ("deep", "deep", "deep", 0.2), ("narrow", "narr", "narrow", 0.2), ("reversed", "rev", "base", -0.2)):
    z = "P_%s_0" % (tag if tag != "rev" else "base"); rP = R("P_%s_2" % tag, z, V); rE = R("E_%s_2" % tag, z, V); s = st[sk]["D2E1"]; pred = K * (4 - s)
    rows[lab] = dict(R_P=rP, R_E=rE, dR=rE - rP, pred=pred, pass_=bool(abs(rE - rP - pred) <= 0.003))
res["P2"] = dict(rows, pass_=bool(all(v["pass_"] for v in rows.values())))
rPs, rEs = R("P_slow_1", "P_base_0", 0.1), R("E_slow_1", "P_base_0", 0.1); d02 = rows["base"]["dR"]
res["P3"] = dict(R_P=rPs, R_E=rEs, dR=rEs - rPs, pred=d02 / 4, pass_=bool(abs(rEs - rPs - d02 / 4) <= 0.0015))
rh = R("E_half_2", "E_half_0", 0.2); rb = R("E_base_2", "E_base_0", 0.2)
res["P4"] = dict(R_half=rh, R_dt01_ownzero=rb, diff=rh - rb, pass_=bool(abs(rh - rb) < 0.001))
rR = R("R_base_2", "P_base_0", 0.2); s = st["base"]["D2E1"]; pr = K * (4 - s * s / 4)
res["P5"] = dict(R=rR, dR=rR - rows["base"]["R_P"], pred=pr, pass_=bool(abs(rR - rows["base"]["R_P"] - pr) <= 0.004 and abs(rR - 1) < 0.02))
Sp = L("spin"); res["P6"] = dict(ratio={k: v["ratio"] for k, v in Sp["runs"].items()},
                                spin_resolved={k: {s: v["H"][s] / v["PS"][s] for s in ("up", "down")} for k, v in Sp["runs"].items()})
res["P6"]["pass_"] = bool(all(abs(v - 0.993) <= 0.003 for v in res["P6"]["ratio"].values()))
res["C1"] = {k: C[k] for k in C if k.startswith("C1")}
res["C2"] = dict(rel=F("E_base_0") / F("P_base_0") - 1, pass_=bool(abs(F("E_base_0") / F("P_base_0") - 1) < 1e-3))
res["C3"] = dict(RJE1=C["C3_RJE1_vs_RJE"], D2E1=C["C3_D2E1_vs_MIX"], pass_=bool(C["C3_RJE1_vs_RJE"] < 1e-14 and C["C3_D2E1_vs_MIX"] < 1e-14))
mz = Sp["mz"]["up"]; res["C4"] = dict(ratio=mz["H"] / mz["PS"], pass_=bool(abs(mz["H"] / mz["PS"] - 0.9930) < 1e-4))
res["core_h"] = sum(json.loads(p.read_text())["_meta"]["wall_s"] for p in OUT.glob("EX-GC64_*.json")) / 3600
(OUT / "summary.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
