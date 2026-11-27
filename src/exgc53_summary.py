"""EX-GC53 scoring from the committed JSON records (separate from the drivers' printouts). Writes results/EX-GC53/summary.json."""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC53"
L = lambda f: json.loads((OUT / f).read_text())["result"]
A, B, C0, CD, E, CK = L("A.json"), L("B.json"), L("C0.json"), L("CD.json"), L("E.json"), L("checks312.json")
s = {}
a = A["a_uniform"]
s["P1"] = dict(a=a["dt"] / a["readout"], weakN=A["a_weakN"]["dt"] / A["a_weakN"]["readout"], b_rel=abs(A["b_gauge"]["dt"]) / a["dt"], c_over_quarter=A["c_one_side"]["dt"] / (a["dt"] / 4),
               pass_a=abs(a["dt"] / a["readout"] - 1) <= 0.02, pass_weakN=abs(A["a_weakN"]["dt"] / A["a_weakN"]["readout"] - 1) <= 0.02,
               pass_b=abs(A["b_gauge"]["dt"]) / a["dt"] < 1e-3, pass_c=abs(A["c_one_side"]["dt"] / (a["dt"] / 4) - 1) <= 0.02)
s["P2"] = dict(rel_total=B["rel_total"], dHdb=B["dHdb_err_max"], pass_=B["rel_total"] <= 1e-10 and B["dHdb_err_max"] <= 1e-8)
lo = CD["loops"]; enc = {h: lo["enclosing"][h] for h in lo["enclosing"]}
p3a = dict(ell0_rel=max(abs(v["0"]["readout"]) / abs(v["2"]["readout"]) for v in enc.values()), flip=max(abs(v["2"]["readout"] + v["-2"]["readout"]) / abs(v["2"]["readout"]) for v in enc.values()),
           N_diff=CD["sym"]["N_diff"])
ds = np.array([16, 24, 32, 40.]); ip = np.array([lo["inplane"][str(int(d))]["2"]["readout"] for d in ds]); ax = np.array([lo["axis"][str(int(d))]["2"]["readout"] for d in ds])
slope_ip = float(np.polyfit(np.log(ds), np.log(np.abs(ip)), 1)[0]); slope_ax = float(np.polyfit(np.log(ds), np.log(np.abs(ax)), 1)[0])
ratios = (ax / ip).tolist()
d_rows = {f"inplane{int(d)}": float(lo["inplane"][str(int(d))]["2"]["readout"] / lo["inplane"][str(int(d))]["dipole"]) for d in ds if d >= 24}
d_rows.update({f"axis{int(d)}": float(lo["axis"][str(int(d))]["2"]["readout"] / lo["axis"][str(int(d))]["dipole"]) for d in ds if d >= 24})
d_rows.update({f"enc{h}": float(v["2"]["readout"] / v["dipole"]) for h, v in enc.items()})
tl = CD["tilted"]
s["P3"] = dict(a=p3a, pass_a=p3a["ell0_rel"] <= 1e-12 and p3a["flip"] <= 1e-6 and p3a["N_diff"] <= 1e-10,
               slope_inplane=slope_ip, slope_axis=slope_ax, pass_b=abs(slope_ip + 3) <= 0.15 and abs(slope_ax + 3) <= 0.15,
               axis_over_inplane=ratios, pass_c=all(abs(r + 2) <= 0.1 for r in ratios),
               over_dipole=d_rows, pass_d=all(abs(v - 1) <= 0.05 for v in d_rows.values()),
               tilted_ratio=tl["ratio_tilted"] / tl["ratio_aligned"], pass_e=abs(tl["ratio_tilted"] / tl["ratio_aligned"] - 1) <= 0.03,
               wave_over_readout={h: v["2"]["wave_energy"] / v["2"]["readout"] for h, v in enc.items()}, lattice_factor=enc["20"]["2"]["lattice_factor"])
R = {k: v["R_dyn"] for k, v in E.items() if "R_dyn" in v}
s["P4"] = dict(R=R, pass_R=all(abs(R[k] - 1) <= 0.02 for k in ("0.1_up", "0.2_up", "0.3_up")), spin_diff=abs(R["0.2_up"] - R["0.2_down"]), pass_spin=abs(R["0.2_up"] - R["0.2_down"]) <= 0.01)
R52 = json.loads((OUT.parent / "EX-GC52" / "D.json").read_text())["result"]["0.2"]["R_dyn"]
lam_star = 4 * (R52 - 1) / (R52 - R["0.2_up"]); s["P5"] = dict(R52=R52, R53=R["0.2_up"], lambda_star=lam_star, pass_=abs(lam_star - 4) <= 0.3)
s["checks"] = dict(C0=C0, pass_C0=C0["ratio_min"] >= 0.97 and C0["ratio_max"] <= 1.03, C1=CK["1"], C2=CK["2"], C3=CK["3"])
ph = OUT / "posthoc_freespace.json"
if ph.exists():
    s["posthoc_freespace"] = json.loads(ph.read_text())
(OUT / "summary.json").write_text(json.dumps(s, indent=1, default=float))
for k, v in s.items():
    print(k, json.dumps(v, default=float)[:700])
