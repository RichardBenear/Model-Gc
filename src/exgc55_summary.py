"""EX-GC55 scoring: results/EX-GC55/*.json -> summary.json. Registered P1-P5, C1-C3 (spec/EX-GC55_plan.md)."""
import json
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1] / "results" / "EX-GC55"
g = lambda n: json.loads((R / f"EX-GC55_{n}.json").read_text())
gam = 1 / np.sqrt(0.96); out = {}

A = g("A"); LF = A["lattice_factor"]
out["P1"] = dict(rel={k: v["rel"] for k, v in A["gauge"].items()}, max_rel=max(v["rel"] for v in A["gauge"].values()))
out["P1"]["pass_"] = out["P1"]["max_rel"] < 1e-3
out["P1"]["diag_cov_oldinit"] = {k: v["diag_cov_oldinit"] for k, v in A["gauge"].items()}; out["P1"]["diag_adv_newinit"] = {k: v["diag_adv_newinit"] for k, v in A["gauge"].items()}
ring = A["ring"]
out["P2"] = dict(uniform=A["a_uniform"]["ratio"] / LF, one_side=A["c_one_side"]["ratio_to_a"] / 0.25, ring={k: v["ratio_over_LF"] for k, v in ring.items()})
out["P2"]["pass_uniform"] = abs(out["P2"]["uniform"] - 1) < 0.005; out["P2"]["pass_one_side"] = abs(out["P2"]["one_side"] - 1) < 0.02
out["P2"]["pass_ring"] = all(abs(v - 1) < 0.01 for v in out["P2"]["ring"].values())

B = {r: g(f"B_{r}") for r in ("P", "MIX", "RJE")}
out["P3"] = {r: dict(drift_rel=B[r]["drift_rel"], Ei0=B[r]["Ei0"]) for r in B}
out["P3"]["pass_"] = B["P"]["drift_rel"] < 1e-4 and B["RJE"]["drift_rel"] < 1e-4 and B["MIX"]["drift_rel"] > 1e-3

C = {f"{r}_{V}": g(f"C_{r}_{V}") for r in ("P", "RJE") for V in ("0.0", "0.2")}
Rr = {r: C[f"{r}_0.2"]["F1"] / C[f"{r}_0.0"]["F1"] * gam for r in ("P", "RJE")}
insitu = {r: {k: C[f"{r}_0.2"]["runs"][k]["j0_tot"][0] / C["P_0.2"]["runs"][k]["j0_tot"][0] for k in ("pair", "one")} for r in ("RJE",)}
out["P4"] = dict(R_P=Rr["P"], R_RJE=Rr["RJE"], diff=Rr["RJE"] - Rr["P"], pass_=abs(Rr["RJE"] - Rr["P"] - 0.030) <= 0.012,
                 RJE_lorentz_2pct=abs(Rr["RJE"] - 1) < 0.02, source_ratio_insitu=insitu, F1={k: v["F1"] for k, v in C.items()})
out["P5"] = dict(R_P=Rr["P"], ratio_to_53=Rr["P"] / 0.9974792, pass_=abs(Rr["P"] / 0.9974792 - 1) < 0.01)

ck = g("checks"); h = g("B_P_half")
out["C1"] = ck["C1"]; out["C2"] = dict(drift_dt=B["P"]["drift_rel"], drift_half=h["drift_rel"], factor=B["P"]["drift_rel"] / h["drift_rel"])
out["C2"]["pass_"] = out["C2"]["factor"] >= 10; out["C3"] = ck["C3"]; out["C3"]["pass_"] = ck["C3"]["herm"] < 1e-12 and ck["C3"]["expect"] < 1e-12
names = ["A", "checks", "B_P", "B_MIX", "B_RJE", "B_P_half"] + [f"C_{k}" for k in C]
out["core_h"] = sum(g(n)["_meta"]["wall_s"] for n in names) / 3600
(R / "summary.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))
