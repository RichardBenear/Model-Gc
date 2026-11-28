"""EX-GC61 post-hoc diagnostics (NOT scored) -> results/EX-GC61/posthoc.json.
(1) D1 with the PREPARATION transformed too: in the shear description the probe starts from U psi0 (the same physical preparation as the
    free psi0), as the Codex 60 review asked; the registered D1 started every run from psi0. Runs: H_g, P shear, sum (sigma = 3), H_g at 2s.
(2) Arrival: continuum retarded point magnetic dipole m(t) = m r(t - d) for the curl of beta: on axis Omega ~ r + d r', in the plane
    Omega ~ r + d r' + d^2 r''; predicted half-rise times against 61B's measured ones; and the first time |Omega| exceeds 1e-3 of settled.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc61 as Q, exgc60 as P60


def d1(kind):
    one = np.ones((48,) * 3); s = 2 * P60.S0 if kind == "Hg2" else P60.S0; Hg, expG = Q.gauge_ops(s, 48); psi = expG(P60.probe(3.0, "x"), +1)
    if kind in ("Hg", "Hg2"):
        r = Q.run_static(one, None, "PS", psi, extra=Hg, tr=expG)
    else:
        r = Q.run_static(one, P60.field("shear" if kind == "P" else "sum", "xy"), "PS", psi, tr=expG)
    (Q.PART / f"ph_{kind}.json").write_text(json.dumps(r))


def ramp_d(u, k):
    tau = Q.TAU; w = np.pi / (2 * tau); u = np.asarray(u)
    if k == 0:
        v = np.sin(w * u) ** 2
    elif k == 1:
        v = w * np.sin(2 * w * u)
    else:
        v = 2 * w * w * np.cos(2 * w * u)
    out = np.where(u < 0, 0.0, np.where(u < tau, v, 1.0 if k == 0 else 0.0)); return out


def arrival():
    B = json.loads((Q.OUT / "EX-GC61_B97.json").read_text()); u = np.linspace(-2, 20, 22001); res = {}
    for key, th in B["t_half"].items():
        kind, d = key.split("_"); d = float(d); f = ramp_d(u, 0) + d * ramp_d(u, 1) + (d * d * ramp_d(u, 2) if kind == "plane" else 0)
        uh = float(u[np.argmax(f >= 0.5)]); om = np.array(B["om"][key]); t = Q.DT * np.arange(len(om)); first = float(t[np.argmax(np.abs(om / om[500]) > 1e-3)])
        res[key] = dict(meas=th, dipole_pred=d + uh, diff=th - (d + uh), first_1e3=first, first_minus_d=first - d)
    return res


if __name__ == "__main__":
    job = sys.argv[1]
    if job != "score":
        d1(job)
    else:
        S = json.loads((Q.OUT / "EX-GC61_D1_3.0.json").read_text()); ang = lambda rec, off=0: np.arctan2(np.array(rec)[:, off + 1], np.array(rec)[:, off])
        fr = ang(S["free"]); sl = lambda a: P60.slope(a - fr, Q.DT) / P60.S0; R = {k: json.loads((Q.PART / f"ph_{k}.json").read_text()) for k in ("Hg", "P", "sum", "Hg2")}
        sigP = abs(ang(S["P"]["raw"])[-1] - fr[-1]); k_P = sl(ang(S["P"]["raw"])); k_Hg = sl(ang(S["Hg"]["raw"])); k_rot = sl(ang(S["rot"]))
        o = dict(max_Hg_tr_over_signal=float(np.abs(ang(R["Hg"]["tr"]) - fr).max() / sigP), max_Hg2_over_Hg=float(np.abs(ang(R["Hg2"]["tr"]) - fr).max() / np.abs(ang(R["Hg"]["tr"]) - fr).max()),
                 k_Ptr=sl(ang(R["P"]["tr"])), k_P_minus_k_Hg=k_P - k_Hg, k_sumtr=sl(ang(R["sum"]["tr"])), k_rot=k_rot, k_Hg_tr=sl(ang(R["Hg"]["tr"])))
        o["c_rel"] = (o["k_Ptr"] - o["k_P_minus_k_Hg"]) / o["k_P_minus_k_Hg"]; o["d_rel"] = (o["k_sumtr"] - (k_rot + o["k_P_minus_k_Hg"])) / abs(k_rot)
        out = dict(D1_matched_prep=o, arrival=arrival()); print(json.dumps(out, indent=1)); (Q.OUT / "posthoc.json").write_text(json.dumps(out, indent=1, default=float))
