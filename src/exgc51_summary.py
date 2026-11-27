"""EX-GC51 scoring and post-run diagnostics, recomputed from the committed JSON records only (separate path from the
drivers' own printouts; shares their JSON outputs). Writes results/EX-GC51/summary.json. EXPLORATORY."""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC51"
L = lambda f: json.loads((OUT / f).read_text())["result"]
A, B, E, D, K = L("A.json"), L("B.json"), L("E.json"), L("D.json"), L("checks123.json")
C = {}; [C.update(L(f)) for f in ("C_m0.25.json", "C_m0.5.json")]
U = sorted(B["rows"], key=float); s = {}
# P1
s["P1"] = dict(c_N=A["field"]["c"], c_light=A["light"]["c"], ratio=A["c_ratio"], pass_cN=abs(A["field"]["c"] - 1) <= 0.02, pass_ratio=abs(A["c_ratio"] - 1) <= 0.01)
# P2 (registered: Geff within 1%; beta_sep as coded) + post-hoc beta separation with strength s = kappa*Q
ge = [v for u in U for v in B["rows"][u]["Geff"].values()]
ph = {}
for d in ("axis", "face", "body"):
    sQ = np.array([B["rows"][u]["kappa"] * B["rows"][u]["info"]["Q"] for u in U]); a1 = np.array([B["rows"][u]["fits"][d]["coef"][0] for u in U]); a2 = np.array([B["rows"][u]["fits"][d]["coef"][1] for u in U])
    cc = np.linalg.lstsq(np.stack([sQ, sQ ** 2, sQ ** 3], 1), a2, rcond=None)[0]; k1 = np.mean(a1 / sQ); ph[d] = float(2 * cc[1] / k1 ** 2)
s["P2"] = dict(Geff_min=min(ge), Geff_max=max(ge), pass_Geff=all(abs(g - 1) <= 0.01 for g in ge), beta_sep=B["beta_sep"],
               pass_beta=all(abs(v["beta"] - 1) <= 0.05 for v in B["beta_sep"].values()),
               beta_raw={u: {d: B["rows"][u]["fits"][d]["beta"] for d in ("axis", "face", "body")} for u in U}, posthoc_beta_sep_kappaQ=ph)
# P3 (A1 scored version)
rows = {k: dict(a_over_semi=v["a_over_semi"], a_over_GR=v["a_over_GR"], semi_over_GR=v["a_semi"] / v["a_GR"], one_plus_rm=v["pred_inf"],
                diag_a_over_semi_full=v["a"] / v["diag_a_semi_full"], diag_plain_over_semi_full=v["diag_a_plain"] / v["diag_a_semi_full"]) for k, v in C.items()}
ratios = {}
for r in (1.0, 0.5, 0.0):
    a, b = C[f"0.25_{r}"], C[f"0.5_{r}"]
    ratios[str(r)] = dict(a_ratio=a["a"] / b["a"], semi_ratio=a["a_semi"] / b["a_semi"], meas_over_semi=(a["a"] / b["a"]) / (a["a_semi"] / b["a_semi"]),
                          GR_ratio=a["a_GR"] / b["a_GR"], diag_plain_over_full=(a["diag_a_plain"] / b["diag_a_plain"]) / (a["diag_a_semi_full"] / b["diag_a_semi_full"]))
s["P3"] = dict(rows=rows, ratios=ratios, pass_a_semi=all(abs(v["a_over_semi"] - 1) <= 0.03 for v in rows.values()),
               pass_ratio=all(abs(v["meas_over_semi"] - 1) <= 0.01 for v in ratios.values()), pass_r0_equal=abs(ratios["0.0"]["a_ratio"] - 1) <= 0.01)
# P4, P5
s["P4"] = dict(momentum_balance=D["momentum_balance"], inertia_balance_diag=D["inertia_balance"], pass_=D["momentum_balance"] < 0.05)
s["P5"] = dict(rows={b: dict(wave=v["ratio_wave"], eikonal=v["ratio_eikonal"]) for b, v in E.items()},
               pass_wave=all(abs(v["ratio_wave"] - 2) <= 0.06 for v in E.values()), pass_eikonal=all(abs(v["ratio_eikonal"] - 2) <= 0.06 for v in E.values()))
s["checks"] = dict(C1=K["C1"]["G_static_vs_49"], C2=K["C2"]["maxres"], C3=K["C3"]["rel"])
(OUT / "summary.json").write_text(json.dumps(s, indent=1, default=float))
for k, v in s.items(): print(k, {kk: vv for kk, vv in v.items() if not isinstance(vv, dict)} if isinstance(v, dict) else v)
print("posthoc beta", ph)
