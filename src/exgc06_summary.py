"""EX-GC06 summary figure: growth rate and frequency of the clocked core's l = 0, 1 modes vs tau_V."""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC06")
runs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "kappa*.json")))]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
sty = {0: "-", 1: "--", 2: ":", 3: "-."}
col = {0.5: "#9bbb59", 0.75: "#8064a2", 0.9: "#3b6ea5", 1.0: "#c0504d"}
summary = []
for r in runs:
    if not r.get("particle"):
        summary.append(dict(kappa=r["kappa"], dv=r["dv"], particle=False))
        continue
    rec = dict(kappa=r["kappa"], dv=r["dv"], particle=True, M=r["M"], W0=r["W0"], lead={})
    for l in range(4):
        pts = sorted([(x["tauV"], x["lead"]) for x in r["rows"] if x["l"] == l and x["lead"]], key=lambda t: t[0])
        rec["lead"][l] = [dict(tauV=t, re=m["re"], im=m["im"]) for t, m in pts]
        if r["dv"] != 0.1 or l > 1:
            continue
        c = col.get(r["kappa"], "k")
        ax[0].plot([t for t, m in pts], [m["re"] for t, m in pts], sty[l], marker="o", ms=3, color=c,
                   label=f"κ = {r['kappa']:g}, l = {l}")
        osc = [(t, m) for t, m in pts if m["im"] > 1e-3]
        ax[1].plot([t for t, m in osc], [m["im"] for t, m in osc], sty[l], marker="o", ms=3, color=c)
    summary.append(rec)
ax[0].axhline(0, color="k", lw=0.6)
for a in ax:
    a.set_xscale("log")
    a.set_xticks([2, 3, 5, 10, 20, 50])
    a.set_xticklabels(["2", "3", "5", "10", "20", "50"])
    a.minorticks_off()
ax[0].set_xlabel("inhibitor delay τ_V"); ax[0].set_ylabel("growth rate of leading core mode")
ax[0].set_title("l = 0 (breathing, solid) goes unstable first; l = 1 (dashed) later", fontsize=9)
ax[0].legend(fontsize=7, ncol=2)
ax[1].set_xlabel("inhibitor delay τ_V"); ax[1].set_ylabel("angular frequency ω (oscillatory modes)")
ax[1].set_title("Clock frequency: ω ≈ 0.3–0.4 at onset (period ≈ 18–20)", fontsize=9)
fig.suptitle("EX-GC06 (EXPLORATORY): a slow core inhibitor gives a breathing (clocked) Gc particle", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC06_clock.png"), dpi=150)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
for s in summary:
    print(s["kappa"], s["dv"], s.get("M"), {l: [(p["tauV"], round(p["re"], 3), round(p["im"], 3)) for p in v]
                                            for l, v in s.get("lead", {}).items() if l < 2})
