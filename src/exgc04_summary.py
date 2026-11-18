"""EX-GC04 summary and figure: bound-pair free fall vs single particle (reads results/EX-GC04/*.json)."""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC04")
L = lambda f: json.load(open(os.path.join(OUT, f)))

s = {h: L(f"single_gradient_h{h}.json") for h in ("0.1", "0.125")}
mu1 = {h: s[h]["mu"] for h in s}
M1 = {h: s[h]["runs"]["0.0001"]["Mc"] for h in s}      # E = +1e-4 single consumption (equals E = -1e-4 to 1e-12)
low = L("single_lowGinf_h0.1.json")
dMdG = (M1["0.1"] - low["runs"]["0.0001"]["Mc"]) / (1 - low["Ginf"])
pairs = []
for f in sorted(glob.glob(os.path.join(OUT, "pair_gradient_d*_h*.json"))):
    p = json.load(open(f))
    h = os.path.basename(f).split("_h")[1][:-5]
    pairs.append(dict(d=p["d"], h=h, fall_ratio=p["mu_pair"] / mu1[h], mass_ratio=p["M_pair"] / (2 * M1[h]),
                      v_rel_E0=p["runs"]["0.0"]["v"][1] - p["runs"]["0.0"]["v"][0], tidal=p["tidal"]))
far = [p for p in pairs if p["h"] == "0.1" and p["d"] >= 14]
A = np.array([[1 / p["d"], 1 / p["d"] ** 2] for p in far])
coef = np.linalg.lstsq(A, np.array([1 - p["fall_ratio"] for p in far]), rcond=None)[0]
vd, excluded = [], []
for f in sorted(glob.glob(os.path.join(OUT, "vd_*.json"))):
    for r in json.load(open(f)):
        r = dict(d=r["d"], v=r["v_right"], ok=r["ok"], src=os.path.basename(f))
        # d = 8.0 from the profile guess converged to an outlier (v = 6.9e-3, 10x its neighbours); continuation from
        # both sides fails for 7.625 <= d <= 8.2, so that window is reported as unresolved, not plotted
        (vd if r["ok"] and not (7.62 < r["d"] < 8.25) else excluded).append(r)
vd.sort(key=lambda r: r["d"])
summary = dict(mu_single=mu1, M_single=M1, dM_dGinf=dMdG, mu_single_lowGinf=low["mu"], Ginf_low=low["Ginf"],
               far_field_prediction_coef=dMdG / (8 * np.pi), fit_a_over_d_plus_b_over_d2=list(coef),
               pairs=pairs, vd=vd, vd_excluded=excluded)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
for p in pairs:
    print("d %.3f h %s fall %.4f mass %.4f pred deficit %.4f" % (p["d"], p["h"], p["fall_ratio"], p["mass_ratio"],
                                                                dMdG / (8 * np.pi * p["d"])))
print("dM/dG", dMdG, "pred coef", dMdG / (8 * np.pi), "fit", coef, "mu low G", low["mu"], "mu", mu1)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for lo, hi in ((0, 7.62), (8.25, 99)):
    seg = [r for r in vd if lo < r["d"] < hi]
    ax[0].plot([r["d"] for r in seg], [r["v"] for r in seg], ".-", color="#3b6ea5", ms=4)
ax[0].axvspan(7.62, 8.25, color="0.85", lw=0)
ax[0].text(7.94, 3e-5, "no steady\npair\nfound", ha="center", fontsize=7)
ax[0].axvline(10.11, color="#c0504d", ls=":", lw=1)
ax[0].axhline(0, color="k", lw=0.6)
ax[0].axvline(7.145, color="#c0504d", ls="--", lw=1)
ax[0].set_yscale("symlog", linthresh=1e-5)
ax[0].set_xlabel("separation d"); ax[0].set_ylabel("speed of each particle (> 0 apart, < 0 together)")
ax[0].set_title("Pair speed v(d): bound at d = 7.145 and 10.11", fontsize=10)
h1 = [p for p in pairs if p["h"] == "0.1"]
dd = np.array([p["d"] for p in h1])
ax[1].plot(dd, [p["fall_ratio"] for p in h1], "o", color="#3b6ea5", label="fall speed, pair / single")
ax[1].plot(dd, [p["mass_ratio"] for p in h1], "s", color="#9bbb59", label="consumption, pair / 2 singles")
x = np.linspace(6.5, 22, 200)
ax[1].plot(x, 1 - dMdG / (8 * np.pi * x), "--", color="#3b6ea5", lw=1, label="far-field estimate 1 − (dM/dG)/(8πd)")
hc = [p for p in pairs if p["h"] == "0.125"]
ax[1].plot([p["d"] for p in hc], [p["fall_ratio"] for p in hc], "x", color="k", label="fall, grid h = 0.125")
ax[1].axvline(7.145, color="#c0504d", ls="--", lw=1)
ax[1].axhline(1, color="k", lw=0.6)
ax[1].set_xlabel("separation d"); ax[1].set_title("Pair in an external G gradient")
ax[1].legend(fontsize=8, loc="lower right")
fig.suptitle("EX-GC04 (EXPLORATORY): a bound pair of identical Model Gc particles falls 6% slower", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC04_boundpair.png"), dpi=150)
