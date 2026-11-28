"""EX-GC54 scoring: results/EX-GC54/*.json (+ 52D, 53E records) -> summary.json. Registered predictions P1-P6, checks C1-C3."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R54 = ROOT / "results" / "EX-GC54"
ld = lambda p: json.loads(Path(p).read_text())
g = lambda name: ld(R54 / f"EX-GC54_{name}.json")
gam = 1 / np.sqrt(1 - 0.04)

E53 = ld(ROOT / "results/EX-GC53/E.json")["result"]; D52 = ld(ROOT / "results/EX-GC52/D.json")["result"]
F53_0, F53_2 = E53["0.0_up"]["F1"], E53["0.2_up"]["F1"]; F52_0, F52_2 = D52["0.0"]["F1"], D52["0.2"]["F1"]
R52 = F52_2 / F52_0 * gam; R53 = F53_2 / F53_0 * gam; slope = (R53 - R52) / 4.0
ops = g("ops"); lam_eff = {"D2q": ops["0.2"]["D2q"]["lam_eff"], "D2E": ops["0.2"]["D2E"]["lam_eff"]}
out = dict(R52=R52, R53=R53, slope=slope, lam_eff=lam_eff)

# P2 operators
out["P2"] = dict(sum_rel={r: ops["0.2"][r]["rel_sum"] for r in lam_eff}, k_rel={r: ops["0.2"][r]["rel_k"] for r in lam_eff})
out["P2"]["pass_sum"] = all(v < 1e-10 for v in out["P2"]["sum_rel"].values()); out["P2"]["pass_k"] = {r: abs(v) < 0.01 for r, v in out["P2"]["k_rel"].items()}

# 54B rules
rows = {}
for r in ("D2q", "D2E"):
    F0, F2 = g(f"B_{r}_0.0")["F1"], g(f"B_{r}_0.2")["F1"]; R = F2 / F0 * gam
    rows[r] = dict(F1_0=F0, F1_02=F2, R_dyn=R, R_linear_pred=R53 + slope * (lam_eff[r] - 4))
out["B"] = rows
out["P3"] = dict(R=rows["D2q"]["R_dyn"], registered=0.977, band=0.008, number_pass=abs(rows["D2q"]["R_dyn"] - 0.977) <= 0.008,
                 lorentz_clause_fails_as_predicted=abs(rows["D2q"]["R_dyn"] - 1) >= 0.02)
out["P4"] = dict(R=rows["D2E"]["R_dyn"], number_pass=abs(rows["D2E"]["R_dyn"] - 1.002) <= 0.008, lorentz_pass=abs(rows["D2E"]["R_dyn"] - 1) < 0.02)

# 54C lambda = 2
C = g("C_P_0.2"); F0i = 0.5 * (F52_0 + F53_0); Rc = C["F1"] / F0i * gam
out["P5"] = dict(F1=C["F1"], F0_interp=F0i, R=Rc, pass_=abs(Rc - 1.056) <= 0.006)
pts = [(0.0, R52), (2.0, Rc), (4.0, R53), (lam_eff["D2E"], rows["D2E"]["R_dyn"]), (lam_eff["D2q"], rows["D2q"]["R_dyn"])]
x, y = np.array(pts).T; co = np.polyfit(x, y, 1); co2 = np.polyfit(x, y, 2)
out["linearity"] = dict(points=pts, fit=co.tolist(), resid=(y - np.polyval(co, x)).tolist(), quad=co2.tolist(), lam_star_lin=float((1 - co[1]) / co[0]))

# P1 no-go
A = g("A"); out["P1"] = dict(max_ratio=A["max_ratio"], pass_=A["max_ratio"] < 0.02)

# P6 smooth gauge
D = g("D"); rel = {k: v["rel"] for k, v in D.items() if isinstance(v, dict) and "Lam" in v}
sc1 = D["64_0.01"]["dt"] / D["64_0.001"]["dt"]; sc2 = D["64_0.001"]["dt"] / D["64_0.0001"]["dt"]
out["P6"] = dict(rel=rel, scale_10x=[sc1, sc2], pass_small=all(D[f"{L}_0.001"]["rel"] < 1e-3 for L in (32, 64, 128)),
                 pass_linear=all(abs(s / 10 - 1) < 0.1 for s in (sc1, sc2)))

# C1: 53's law through the 54 code vs 53E stored P1(t), t <= 8
C1 = g("C1_P_0.2"); t = np.array(C1["runs"]["pair"]["t"]); P1 = np.array([p[0] for p in C1["runs"]["pair"]["Py"]]) - np.array([p[0] for p in C1["runs"]["one"]["Py"]])
t53 = np.array(E53["0.2_up"]["t"]); P53 = np.array(E53["0.2_up"]["P1"]); m = t53 <= 8 + 1e-9
out["C1"] = dict(maxdiff_rel=float(np.max(np.abs(P1 - P53[m][:len(P1)])) / np.max(np.abs(P53[m]))))
out["C1"]["pass_"] = out["C1"]["maxdiff_rel"] < 1e-10

# C2: D2-E dt/2 vs dt, 2 <= t <= 8
C2 = g("C2_D2E_0.2"); B2 = g("B_D2E_0.2")
def P1of(rec):
    t = np.array(rec["runs"]["pair"]["t"]); return t, np.array([p[0] for p in rec["runs"]["pair"]["Py"]]) - np.array([p[0] for p in rec["runs"]["one"]["Py"]])
th, Ph = P1of(C2); tf, Pf = P1of(B2); Phi = np.interp(tf, th, Ph); m = (tf >= 2) & (tf <= 8 + 1e-9)
out["C2"] = dict(maxdiff_rel=float(np.max(np.abs(Phi[m] - Pf[m])) / np.max(np.abs(Pf[m])))); out["C2"]["pass_"] = out["C2"]["maxdiff_rel"] < 1e-5

# C3 ring
ring = g("ring"); out["ring"] = ring
c3 = {}
for r in ("D2q", "D2E"):
    Jp, Jm = np.array(ring["2"][r]["J"]), np.array(ring["-2"][r]["J"])
    c3[r] = dict(ell0_max=ring["0"][r]["maxabs"], flip=float(np.abs(Jp + Jm).max() / np.abs(Jp).max()), Jz=float(Jp[2]), Jz_over_P=float(Jp[2] / ring["2"]["J_P"][2]))
out["C3"] = c3; out["C3_pass"] = all(v["ell0_max"] < 1e-15 and v["flip"] < 1e-10 for v in c3.values())

core_s = sum(g(n)["_meta"]["wall_s"] for n in ("B_D2q_0.0", "B_D2q_0.2", "B_D2E_0.0", "B_D2E_0.2", "C_P_0.2", "C1_P_0.2", "C2_D2E_0.2", "A", "D", "ring", "ops"))
PAUSE_S = 1960.0   # C_P_0.2 and C2_D2E_0.2 were paused (SIGSTOP) at the 2 core-h cap from about 23:04:50 to 23:37:30 UTC; wall includes it
out["core_h_wall_sum"] = core_s / 3600; out["core_h"] = (core_s - 2 * PAUSE_S) / 3600
# post hoc (not scored): in-situ source totals at t = 0 (in the packets' own N wells) and the linear prediction with them
P0 = g("C1_P_0.2")["runs"]
for r in ("D2q", "D2E"):
    j = g(f"B_{r}_0.2")["runs"]; le = {k: 4 * j[k]["j0_tot"][0] / P0[k]["j0_tot"][0] for k in j}
    out["B"][r]["posthoc_lam_insitu"] = le; out["B"][r]["posthoc_R_pred_insitu"] = R53 + slope * (le["pair"] - 4)
(R54 / "summary.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))
