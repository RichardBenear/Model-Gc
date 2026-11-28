"""EX-GC57 scoring -> results/EX-GC57/summary.json (registered P1-P6, C1-C4 of spec/EX-GC57_plan.md)."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results" / "EX-GC57"; g = lambda n: json.loads((R / f"EX-GC57_{n}.json").read_text())
A = g("A"); B = g("B"); out = {}
rows = A["rows"]; pa = [r["a"]["E_over_eps"] for r in rows]; pe = [r["a"]["env_over_u2"] for r in rows]
b16 = {r["k0"]: r["b"]["E_over_56pred"] for r in rows if r["width"] == 16.0}
out["P1"] = dict(a_E_over_eps=pa, a_env_over_u2=pe, b_w16=b16, ref56={0.25: 0.9942, 0.5: 0.9972})
out["P1"]["pass_"] = max(abs(x - 1) for x in pa + pe) < 1e-3 and abs(b16[0.25] - 0.9942) < 5e-4 and abs(b16[0.5] - 0.9972) < 5e-4
P = B["positions"]; keys = ["axis_24", "plane_24", "axis_32", "plane_32"]; p2 = {}
for k in keys:
    o = P[k]["2"]; om = P[k]["-2"]; cz = o["Omega_curl"][2]
    odd = lambda a, b: abs(a - b) / 2 / max(abs(a + b) / 2, 1e-300)
    p2[k] = dict(orb_over_curl=o["Omega_orb"] / cz, light_over_curl=o["light"]["Omega_light"] / cz, orb_over_light=o["Omega_orb"] / o["light"]["Omega_light"],
                 flip_orb=odd(o["Omega_orb"], om["Omega_orb"]), flip_light=odd(o["light"]["dt"], om["light"]["dt"]))
out["P2"] = dict(cases=p2)
out["P2"]["pass_orb"] = all(abs(v["orb_over_curl"] - 1) < 0.1 for v in p2.values())
out["P2"]["pass_light"] = all(abs(v["light_over_curl"] - 1) < 0.1 for v in p2.values())
out["P2"]["pass_orb_light"] = all(abs(v["orb_over_light"] - 1) < 0.1 for v in p2.values())
out["P2"]["pass_flip"] = all(min(v["flip_orb"], v["flip_light"]) > 1e6 for v in p2.values())
out["P2"]["pass_"] = all(out["P2"][k] for k in ("pass_orb", "pass_light", "pass_orb_light", "pass_flip"))
r3 = P["axis_32"]["2"]["Omega_orb"] / P["plane_32"]["2"]["Omega_orb"]; rc = P["axis_32"]["2"]["Omega_curl"][2] / P["plane_32"]["2"]["Omega_curl"][2]
out["P3"] = dict(orb_ratio=r3, curl_ratio=rc, pass_=abs(r3 / -2 - 1) < 0.15)
out["P4"] = {k: P[k]["2"]["Omega_torque"] / P[k]["2"]["Omega_orb"] for k in ("axis_24", "axis_32")}; out["P4"]["pass_"] = all(abs(v - 1) < 0.1 for k, v in out["P4"].items() if k != "pass_")
cB = B["cB"]; p5 = {}
for k in ("axis_24", "axis_32"):
    o = P[k]["2"]; cz = o["Omega_curl"][2]
    p5[k] = {s: dict(ratio=o[f"Omega_spin_{s}"] / cz, two_cB=2 * cB[f"{s}.0"][0], rel=o[f"Omega_spin_{s}"] / cz / (2 * cB[f"{s}.0"][0])) for s in ("3", "5")}
out["P5"] = dict(cases=p5, cB=cB)
out["P5"]["pass_"] = all(v[s]["ratio"] < 0.3 and abs(v[s]["rel"] - 1) < 0.25 for v in p5.values() for s in ("3", "5")) and all(v["5"]["ratio"] < v["3"]["ratio"] for v in p5.values())
gr = {}
for k in ("axis_24", "plane_24"):
    gc = P[k]["grad"]; o = P[k]["2"]; gr[k] = dict(orb=abs(gc["Omega_orb"] / o["Omega_orb"]), light=abs(gc["light"]["dt"] / o["light"]["dt"]))
    if k == "axis_24":
        gr[k].update({f"spin_{s}": abs(gc[f"Omega_spin_{s}"] / o[f"Omega_spin_{s}"]) for s in ("3", "5")})
out["P6"] = dict(cases=gr, pass_=all(x < 0.01 for v in gr.values() for x in v.values()))
out["C1"] = dict(A["C1"], pass_=A["C1"]["zero_drift"] < 1e-12 and A["C1"]["rel"] < 1e-4)
out["C2"] = dict(B["C2"], pass_=B["C2"]["rel"] < 1e-12)
out["C3"] = dict(B["C3"], pass_=abs(B["C3"]["uniform_over_2sum"] - 1) < 0.01)
out["C4"] = dict(q_diff=B["C4_q_diff"], N_diff=B["C4_N_diff"], B0_max=B["C4_B0_max"], pass_=B["C4_q_diff"] < 1e-12 and B["C4_N_diff"] < 1e-12 and B["C4_B0_max"] < 1e-15)
out["core_h"] = (A["_meta"]["wall_s"] + B["_meta"]["wall_s"]) / 3600
(R / "summary.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps({k: (v.get("pass_") if isinstance(v, dict) else v) for k, v in out.items()}))
