"""EX-GC61 scoring -> results/EX-GC61/summary.json. Gates exactly as in spec/EX-GC61_plan.md. EXPLORATORY. Our construction; not
LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc61 as Q
import exgc60 as P60

OUT = Q.OUT; L = lambda k: json.loads((OUT / f"EX-GC61_{k}.json").read_text()); s0 = P60.S0; dt = Q.DT
A, B97, B129, C3 = L("A"), L("B97"), L("B129"), L("C3"); res = {}; core = 0.0
for k in ("A", "B97", "B97w", "B129", "C97up", "C97dn", "C129up", "C129dn", "C3", "D1_3.0", "D1_4.0", "D1_5.0", "D2", "D3"):
    core += L(k)["_meta"]["wall_s"]
res["core_h"] = core / 3600 + 109 / 3600 * 0  # C3 counted in its own json

# P1
res["P1"] = dict(diff=A["diff_dt0.1"], ratio=A["ratio"], shuffle=max(A["shuffle_seed1"], A["shuffle_seed2"]))
res["P1"]["pass_"] = bool(A["diff_dt0.1"] < 1e-3 and 3 <= A["ratio"] <= 5 and res["P1"]["shuffle"] < 1e-12)
# P2
p2 = {}
for kind in ("axis", "plane"):
    ds = np.array([8, 12, 16, 24.0]); th = np.array([B97["t_half"][f"{kind}_{int(d)}"] for d in ds]); off = th - ds
    p2[kind] = dict(t_half=th.tolist(), offset=off.tolist(), spread=float(np.abs(off - off.mean()).max()), slope=float(np.polyfit(ds, th, 1)[0]))
res["P2"] = dict(p2, pass_=bool(all(v["spread"] <= 0.75 and abs(v["slope"] - 1) <= 0.03 for v in p2.values())))
# P3
res["P3"] = dict(rms=B97["rms_vs_static193"], rms_vs_zero_box97=B97["rms_vs_static97"], pass_=bool(B97["rms_vs_static193"] < 0.01))
# P4
o97 = np.array(B97["om"]["axis_16"]); o129 = np.array(B129["om"]["axis_16"]); t = dt * np.arange(len(o97)); ref = abs(o129[500]); rel = np.abs(o97 - o129) / ref
early = float(rel[t <= 70 + 1e-9].max()); first = float(t[np.argmax(rel > 1e-3)]) if np.any(rel > 1e-3) else None
res["P4"] = dict(max_rel_t_le_70=early, first_exceed_1e3=first, pass_=bool(early <= 1e-4 and first is not None and 76 <= first <= 92))
# P5-P7
def rates(nb, key):
    u = L(f"C{nb}up")["probes"][key]; d = L(f"C{nb}dn")["probes"][key]; return {"S": Q.odd_rate(u, d, 0), "J": Q.odd_rate(u, d, 6), "L": Q.odd_rate(u, d, 3)}
r40, r90, r129 = rates(97, "axis_16@40"), rates(97, "axis_16@90"), rates(129, "axis_16@40")
gr = C3["Omega_gr"]; m59 = C3["ref"]
res["P5"] = dict(rate=r40, over_GR={k: r40[k] / gr for k in ("S", "J")}, over_59static={k: r40[k] / m59 for k in ("S", "J")},
                 n129_rel={k: r129[k] / r40[k] - 1 for k in ("S", "J")})
res["P5"]["pass_"] = bool(all(abs(res["P5"]["over_GR"][k] - 1) <= 0.015 and 1.010 <= res["P5"]["over_59static"][k] <= 1.035 and abs(res["P5"]["n129_rel"][k]) <= 0.002 for k in ("S", "J")))
Bo = 0.5 * (np.load(Q.PART / "B45_97_up_diag.npy") - np.load(Q.PART / "B45_97_dn_diag.npy")); M = Q.wgrad(Bo); pred, rot_only = Q.pred_xy(M); rd = rates(97, "diag@40")
res["P6"] = dict(rate=rd, M=M.tolist(), pred=pred, ratio={k: rd[k] / pred[k] for k in ("S", "J")}, J_rot_only_ratio=rd["J"] / rot_only["J"],
                 S_rot_only_ratio=rd["S"] / rot_only["S"])
res["P6"]["pass_"] = bool(all(abs(res["P6"]["ratio"][k] - 1) <= 0.01 for k in ("S", "J")) and abs(res["P6"]["J_rot_only_ratio"] - 1) <= 0.01)
# P6 extra (reported): axis prediction from the dynamic field's gradients
Ba = 0.5 * (np.load(Q.PART / "B45_97_up_axis_16.npy") - np.load(Q.PART / "B45_97_dn_axis_16.npy")); pa, _ = Q.pred_xy(Q.wgrad(Ba))
res["axis_formula_ratio"] = {k: r40[k] / pa[k] for k in ("S", "J")}
res["P7"] = dict(rate90=r90, rel_change_J=r90["J"] / r40["J"] - 1, rel_change_S=r90["S"] / r40["S"] - 1, pass_=bool(abs(r90["J"] / r40["J"] - 1) > 0.005))
# P8
ang = lambda rec, off=0: np.arctan2(np.array(rec)[:, off + 1], np.array(rec)[:, off])
def d1(sig):
    D = L(f"D1_{sig}"); fr = ang(D["free"]); slope = lambda a: P60.slope(a - fr, dt) / s0
    o = dict(k_P=slope(ang(D["P"]["raw"])), k_Hg=slope(ang(D["Hg"]["raw"])), k_Ptr=slope(ang(D["P"]["tr"])),
             sigP_end=float(abs(ang(D["P"]["raw"])[-1] - fr[-1])), max_Hg_tr=float(np.abs(ang(D["Hg"]["tr"]) - fr).max()))
    frJ = ang(D["free"], 6); o["J"] = dict(k_P=P60.slope(ang(D["P"]["raw"], 6) - frJ, dt) / s0, k_Hg=P60.slope(ang(D["Hg"]["raw"], 6) - frJ, dt) / s0,
                                         k_Ptr=P60.slope(ang(D["P"]["tr"], 6) - frJ, dt) / s0)
    if sig == "3.0":
        o["k_rot"] = slope(ang(D["rot"])); o["k_sumtr"] = slope(ang(D["sum"]["tr"])); o["k_sumraw"] = slope(ang(D["sum"]["raw"]))
        o["max_Hg2_tr"] = float(np.abs(ang(D["Hg2"]["tr"]) - fr).max())
    o["frac"] = (o["k_P"] - o["k_Hg"]) / o["k_P"]; return o
d3, d4, d5 = d1("3.0"), d1("4.0"), d1("5.0"); diff = d3["k_P"] - d3["k_Hg"]
p8 = dict(s3=d3, s4=d4, s5=d5)
p8["a"] = d3["max_Hg_tr"] / d3["sigP_end"]; p8["b"] = d3["max_Hg2_tr"] / d3["max_Hg_tr"]; p8["c"] = (d3["k_Ptr"] - diff) / diff
p8["d"] = (d3["k_sumtr"] - (d3["k_rot"] + diff)) / abs(d3["k_rot"]); p8["e"] = dict(r4=d4["frac"] / d3["frac"], r5=d5["frac"] / d3["frac"])
p8["pass"] = dict(a=bool(p8["a"] < 1e-3), b=bool(3 <= p8["b"] <= 5), c=bool(abs(p8["c"]) < 0.05), d=bool(abs(p8["d"]) < 1e-3),
                  e=bool(abs(p8["e"]["r4"] / 0.5625 - 1) <= 0.25 and abs(p8["e"]["r5"] / 0.36 - 1) <= 0.25))
p8["pass_"] = bool(all(p8["pass"].values())); res["P8"] = p8
# P9
D2 = L("D2"); A5 = json.loads((Q.ROOT / "results/EX-GC60/EX-GC60_A5.json").read_text())["rows"]["xy_shear_PS"]["S"]
res["P9"] = dict(P_vs_PS=abs(D2["P"]["S"] - D2["PS"]["S"]) / abs(D2["PS"]["S"]), k64=D2["PS"]["S"], k48=A5, rel=D2["PS"]["S"] / A5 - 1)
res["P9"]["pass_"] = bool(res["P9"]["P_vs_PS"] < 1e-7 and abs(res["P9"]["rel"]) < 1e-3)
# P10
D3 = L("D3")
def vec_pred(Mg, v0, k):
    Om = -0.5 * np.array([Mg[2, 1] - Mg[1, 2], Mg[0, 2] - Mg[2, 0], Mg[1, 0] - Mg[0, 1]]); Sy = 0.5 * (Mg + Mg.T); u = v0 / np.linalg.norm(v0)
    p = -K60k(k, "rot") * np.cross(Om, v0) + K60k(k, "sh") * (Sy @ v0 - (u @ Sy @ u) * v0); return p - (p @ u) * u
K60k = lambda k, c: Q.K60[k][c]
def vslope(rec, off):
    r = np.array(rec)[:, off:off + 3]; tt = dt * np.arange(len(r)); return np.array([np.polyfit(tt, r[:, i], 1)[0] for i in range(3)])
p10 = {}
for case in ("flat", "src"):
    for k, off in (("S", 0), ("J", 6)):
        if case == "flat":
            meas = (vslope(D3["flat_G"], off) - vslope(D3["flat_free"], off)) / s0; Mg = Q.GMAT; v0 = np.array(D3["flat_free"][0][off:off + 3])
        else:
            meas = 0.5 * (vslope(D3["src_up"], off) - vslope(D3["src_down"], off)); Mg = np.array(D3["M_src"]); v0 = np.array(D3["src_up"][0][off:off + 3])
        u = v0 / np.linalg.norm(v0); mp = meas - (meas @ u) * u; pp = vec_pred(Mg, v0, k)
        p10[f"{case}_{k}"] = dict(meas_perp=mp.tolist(), pred_perp=pp.tolist(), err=float(np.linalg.norm(mp - pp) / np.linalg.norm(pp)))
res["P10"] = dict(p10, pass_=bool(all(v["err"] <= 0.03 for v in p10.values())))
# checks
res["C1"] = dict(v=A["C1_zero"], pass_=bool(A["C1_zero"] == 0)); res["C2"] = dict(v=B97["C2_energy"]["spread_rel"], pass_=bool(B97["C2_energy"]["spread_rel"] < 1e-4))
res["C3"] = dict(v=C3["rel"], pass_=bool(abs(C3["rel"]) < 1e-10)); res["C4"] = dict(div_rel=B97["C4_div_rel"])
Bw = L("B97w"); res["weighted_variant"] = dict(t_half=Bw["t_half"], rms_vs_own_static=Bw["rms_vs_own_static"])
(OUT / "summary.json").write_text(json.dumps(res, indent=1, default=float))
print({k: v.get("pass_") for k, v in res.items() if isinstance(v, dict) and "pass_" in v}); print("core_h", res["core_h"])
