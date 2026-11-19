"""EX-GC05 summary figure: leading eigenvalues of the Gc1 particle's round and non-round modes vs response times."""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC05")


def lead(modes, l):
    """Leading mode, skipping the l = 1 translation mode and the slow G-diffusion ladder (|re| < 0.02, real)."""
    ms = [m for m in modes if not (abs(m["im"]) < 1e-6 and m["re"] > -0.02 and l == 1)]
    return ms[0] if ms else None


rows = []
for f in sorted(glob.glob(os.path.join(OUT, "scan_*.json"))):
    rows += json.load(open(f))
summary = {}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
cols = {0: "#3b6ea5", 1: "#c0504d", 2: "#9bbb59", 3: "#8064a2"}
for typ, mk in (("ref", "o"), ("GY_k0.2", "s"), ("GY_k0.4", "^")):
    for l in range(4):
        rr = sorted([r for r in rows if r.get("type", "ref") == typ and "tauY" in r and r["l"] == l],
                    key=lambda r: r["tauY"])
        pts = [(r["tauY"], lead(r["core_modes"], l)) for r in rr]
        pts = [(t, m) for t, m in pts if m]
        summary[f"{typ}_l{l}"] = [dict(tauY=t, **m) for t, m in pts]
        if not pts:
            continue
        ax[0].plot([t for t, m in pts], [m["re"] for t, m in pts], mk + "-", color=cols[l], ms=4, lw=0.8,
                   label=f"l = {l}" if typ == "ref" else None, mfc="none" if typ != "ref" else cols[l])
        osc = [(t, m) for t, m in pts if abs(m["im"]) > 1e-3]
        ax[0].plot([t for t, m in osc], [m["re"] for t, m in osc], "*", color="k", ms=8)
ax[0].axvline(140 ** 2 / 81 / 28, color="0.4", ls="--", lw=1)
ax[0].text(8.5, 0.15, "ether Hopf\nτ_Y = 8.64", ha="right", fontsize=8)
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xlabel("Y response time τ_Y"); ax[0].set_ylabel("leading core-mode growth rate")
ax[0].set_title("Filled: G feeds X; open: G feeds Y (k 0.2, 0.4); ★ oscillatory", fontsize=9)
ax[0].legend(fontsize=8)
for l in (0, 2, 3):
    rr = sorted([r for r in rows if "tauW" in r and r["l"] == l], key=lambda r: r["tauW"])
    pts = [(r["tauW"], r["core_modes"][0]) for r in rr if r["core_modes"]]
    summary[f"tauW_l{l}"] = [dict(tauW=t, **m) for t, m in pts]
    ax[1].plot([t for t, m in pts], [m["re"] for t, m in pts], "o-", color=cols[l], label=f"l = {l}")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].axvline(28.3, color="0.4", ls="--", lw=1)
ax[1].text(28.8, -0.07, "drift threshold\n(EX-GC02)", fontsize=8)
ax[1].set_xlabel("core response time τ_W"); ax[1].set_title("Slow core: modes slow down but stay real", fontsize=9)
ax[1].legend(fontsize=8)
fig.suptitle("EX-GC05 (EXPLORATORY): no oscillatory non-round mode, so no rotating Gc1 particle", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC05_rotation_search.png"), dpi=150)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
for k, v in summary.items():
    print(k, [(round(x.get("tauY", x.get("tauW")), 2), round(x["re"], 4), round(x["im"], 3)) for x in v])
