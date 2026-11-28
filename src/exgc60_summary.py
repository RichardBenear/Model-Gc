"""EX-GC60 scoring -> results/EX-GC60/summary.json (registered P1-P8, C1-C3 of spec/EX-GC60_plan.md)."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results" / "EX-GC60"; g = lambda n: json.loads((R / f"EX-GC60_{n}.json").read_text())
A = {s: g(f"A{s}") for s in (3, 4, 5)}; B = g("B"); C = g("C"); D = g("D"); out = {}
rot = [A[s]["rows"][f"{p}_rot_PS"][k] for s in A for p in ("xy", "yz") for k in ("S", "J")]
out["P1"] = dict(rot=rot, links_rot_S=A[3]["rows"]["xy_rot_P"]["S"]); out["P1"]["pass_"] = all(abs(x + 1) < 0.01 for x in rot) and abs(out["P1"]["links_rot_S"] + 0.10) < 0.02
ks = {s: A[s]["rows"]["xy_shear_PS"]["S"] for s in A}; f = {s: A[s]["f"] for s in A}; q = {s: ks[s] / f[s] for s in A}
pp = {s: abs(A[s]["rows"]["xy_shear_P"]["S"] / ks[s] - 1) for s in A}
out["P2"] = dict(k_sh=ks, f=f, k_over_f=q, spread=(max(q.values()) - min(q.values())) / np.mean(list(q.values())), P_vs_PS=pp)
out["P2"]["pass_"] = ks[3] > ks[4] > ks[5] and out["P2"]["spread"] < 0.1 and max(pp.values()) < 1e-6
out["P3"] = {s: abs(A[s]["rows"]["xy_shear_PS"]["J"] / ks[s]) for s in A}; out["P3"]["pass_"] = all(v < 0.1 for k, v in out["P3"].items() if k != "pass_")
add = [abs(A[s]["rows"]["xy_sum_PS"][k] / (A[s]["rows"]["xy_rot_PS"][k] + A[s]["rows"]["xy_shear_PS"][k]) - 1) for s in A for k in ("S", "J")]
out["P4"] = dict(max_dev=max(add), pass_=max(add) < 0.005)
sym = [abs(A[s]["rows"][f"yz_{kd}_PS"][k] / A[s]["rows"][f"xy_{kd}_PS"][k] - 1) for s in A for kd in ("rot", "shear", "sum") for k in ("S", "J")]
out["P5"] = dict(max_dev=max(sym), pass_=max(sym) < 0.01)
out["P6"] = dict(a=B["Hg_vs_gauge_rel"], Hg_over_P=B["Hg_over_P"], pass_a=B["Hg_vs_gauge_rel"] < 1e-3, pass_b=abs(B["Hg_over_P"]) < 0.1); out["P6"]["pass_"] = out["P6"]["pass_a"] and out["P6"]["pass_b"]
out["P7"] = dict(diag=C["diag"]["ratio"], pass_=all(abs(v - 1) < 0.03 for v in C["diag"]["ratio"].values()))
out["P8"] = dict(rel=D["rel"], pass_=abs(D["rel"]) < 0.005)
out["C1"] = dict(v=A[3]["C1"], pass_=A[3]["C1"] == 0.0)
out["C2"] = dict(ratio=A[3]["rows"]["xy_shear2_PS"]["S"] / ks[3], pass_=abs(A[3]["rows"]["xy_shear2_PS"]["S"] / ks[3] / 2 - 1) < 1e-3)  # rates are per unit S0, so 2s gives ratio 2
m59 = json.loads((ROOT / "results/EX-GC59/EX-GC59_main.json").read_text())["runs"]["axis_16"]["Omega_meas"]
out["C3"] = dict(v60=C["axis_16"]["meas"]["S"], v59=m59, rel=abs(C["axis_16"]["meas"]["S"] / m59 - 1)); out["C3"]["pass_"] = out["C3"]["rel"] < 1e-10
out["consistency"] = dict(axis16=C["axis_16"]["ratio"], plane16_pred=C["plane_16"]["pred"])
out["core_h"] = sum(g(n)["_meta"]["wall_s"] for n in ("A3", "A4", "A5", "B", "C", "D")) / 3600
(R / "summary.json").write_text(json.dumps(out, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o)))
print({k: (v.get("pass_") if isinstance(v, dict) and "pass_" in v else "-") for k, v in out.items()}); print(json.dumps(out["P2"], default=float)); print(out["P3"], out["C1"], out["C3"], out["core_h"])
