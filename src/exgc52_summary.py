"""EX-GC52 scoring and post-run diagnostics from the committed JSON records and the saved z = c field slices
(partial_fields_B.npy, SHA-256 committed). Writes results/EX-GC52/summary.json. EXPLORATORY."""
import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC52"
L = lambda f: json.loads((OUT / f).read_text())["result"]
A, BC = L("A.json"), L("BC.json"); s = {}
s["P1"] = dict(max_rel=max(abs(v["rel"]) for v in A.values()), pass_=all(abs(v["rel"]) <= 1e-6 for v in A.values()),
               reported_Q_ratio={k: dict(QV_over_Q0=v["Q"] / A[f"{v['r']}_0.0"]["Q"], gamma_1pV2=(1 + v["V"] ** 2) / np.sqrt(1 - v["V"] ** 2)) for k, v in A.items()})
P2 = {V: dict(amp_rho=BC[V]["amp_rho_mean"], amp_rho_z=BC[V]["amp_rho_z_mean"], contraction_over_gamma=BC[V]["contraction_over_gamma"]) for V in ("0.1", "0.2", "0.3")}
s["P2"] = dict(rows=P2, pass_=all(abs(v["amp_rho"] - 1) <= 0.02 and abs(v["contraction_over_gamma"] - 1) <= 0.02 for v in P2.values()))
P3 = {V: dict(R=BC[V]["R"], pred=BC[V]["R_pred_lattice"], cont=BC[V]["R_cont"], R_over_pred=BC[V]["R"] / BC[V]["R_pred_lattice"], excess_over_4V2=(BC[V]["R"] - 1) / (4 * float(V) ** 2)) for V in ("0.1", "0.2", "0.3")}
s["P3"] = dict(rows=P3, pass_match=all(abs(v["R_over_pred"] - 1) <= 0.01 for v in P3.values()), pass_excess=all(P3[V]["excess_over_4V2"] > 0.5 for V in ("0.2", "0.3")))
s["P5"] = dict(rest_impulse_ratio={V: BC[V]["rest_impulse_ratio"] for V in ("0.1", "0.2", "0.3")}, pass_=abs(BC["0.2"]["rest_impulse_ratio"] - 1) <= 0.03)
# post hoc: the moving spin-up source's charge centroid is displaced in -y (spin-dependent side shift); symmetrise +y/-y
F = np.load(OUT / "partial_fields_B.npy"); n = F.shape[1]; c = n // 2; rr = np.arange(12, 37); ph = {}
for i, V in enumerate(("0.0", "0.1", "0.2", "0.3")):
    psi = np.log(F[i]); r = BC[V]; ref = -r["kappa"] * r["Q"] * r["gamma"] / (4 * np.pi)
    up, dn = np.mean(psi[c, c + rr] * rr) / ref, np.mean(psi[c, c - rr] * rr) / ref
    gp, gm = (F[i][c, c + 17] - F[i][c, c + 15]) / 2, (F[i][c, c - 15] - F[i][c, c - 17]) / 2
    ph[V] = dict(amp_plus_y=up, amp_minus_y=dn, amp_sym=(up + dn) / 2, sym_grad_factor=float((gp - gm) / 2 / gp))
for V in ("0.1", "0.2", "0.3"):
    ph[V]["R_sym"] = BC[V]["R"] * ph[V]["sym_grad_factor"] / ph["0.0"]["sym_grad_factor"]; ph[V]["R_sym_over_pred"] = ph[V]["R_sym"] / BC[V]["R_pred_lattice"]
s["posthoc_symmetrised"] = ph
for f, key in (("D.json", "D"), ("checks12.json", "checks12"), ("checks3.json", "checks3")):
    if (OUT / f).exists():
        s[key] = L(f)
if "D" in s and "0.2" in s["D"]:
    Rd = s["D"]["0.2"]["R_dyn"]; s["P4"] = dict(R_dyn=Rd, R_semi=BC["0.2"]["R"], ratio=Rd / BC["0.2"]["R"], pass_=abs(Rd / BC["0.2"]["R"] - 1) <= 0.03,
                                               posthoc_ratio_to_R_sym=Rd / ph["0.2"]["R_sym"])
    for k in ("t", "P1"):
        [s["D"][v].pop(k, None) for v in s["D"]]
if "checks3" in s and "D" in s:
    s["C3"] = dict(F1_dt=s["D"]["0.2"]["F1"], F1_dt_half=s["checks3"]["C3"]["0.2"]["F1"], rel=s["checks3"]["C3"]["0.2"]["F1"] / s["D"]["0.2"]["F1"] - 1)
    c3 = s["checks3"]["C3"]["0.2"]; d = s["D"]["0.2"]
    s["C3"].update(kappa_dt=d["kappa"], kappa_dt_half=c3["kappa"], posthoc_rel_per_kappa=(c3["F1"] / c3["kappa"]) / (d["F1"] / d["kappa"]) - 1,
                   note="check bug: c3() ran part_D with Vs=(0.2,), which set kappa from the V = 0.2 packet's charge instead of the V = 0 packet's")
    s.pop("checks3")
(OUT / "summary.json").write_text(json.dumps(s, indent=1, default=float))
for k, v in s.items():
    print(k, json.dumps(v, default=float)[:600])
