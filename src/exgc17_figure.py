"""EX-GC17 figure: drag per unit, terminal fall speed (EX-GC16 push / drag) by type, and acceleration response."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
s17 = json.load(open(os.path.join(ROOT, "results", "EX-GC17", "summary.json")))
s16 = json.load(open(os.path.join(ROOT, "results", "EX-GC16", "summary.json")))
kinds = ["single", "double_site", "spread7", "pair_z", "pair_x"]
labels = ["single", "2x rate\none site", "spread\n7 sites", "two units\nstacked", "two units\nside by side"]


def g(k, m0, u, w):
    return next(x for x in s17 if x["kind"] == k and x["m0"] == m0 and x["u"] == u and x["w"] == w)


fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
for j, (m0, col) in enumerate(((50, "#1f77b4"), (500, "#d62728"))):
    u = 0.1 if m0 == 50 else 0.3
    xx = [i + (j - 0.5) * 0.25 for i in range(5)]
    ax[0].plot(xx, [g(k, m0, u, 0)["gamma"] / (m0 / 6) for k in kinds], "o", color=col, label=f"measured, m0 = {m0}")
    ax[0].plot(xx, [g(k, m0, u, 0)["gamma_pred"] / (m0 / 6) for k in kinds], "_", ms=14, mew=2, color=col, alpha=.6,
               label=f"mean field, m0 = {m0}")
    f16 = {x["kind"]: x for x in s16 if x["m0"] == m0 and x["E"] == 0.02}
    ax[1].plot(xx, [-f16[k]["fall_z"] / g(k, m0, u, 0)["gamma"] / 0.04 for k in kinds], "o", color=col,
               label=f"measured, m0 = {m0}")
    ax[1].plot(xx, [-f16[k]["fall_z"] * 6 / m0 / 0.04 for k in kinds], "x", color=col, alpha=.5,
               label=f"push alone (EX-GC16), m0 = {m0}")
for a in ax[:2]:
    a.axhline(1, color="k", lw=.8, ls="--")
    a.set_xticks(range(5), labels, fontsize=8)
    a.legend(fontsize=7)
ax[0].set_ylabel("drag per unit / (m0/6)")
ax[0].set_title("Drag from meeting tokens head-on")
ax[1].set_ylabel("terminal fall speed / 2E")
ax[1].set_title("Terminal fall speed = push / drag (E = 0.02)")
for w, mk in ((0.2, "s"), (1.0, "o")):
    ks = ["single", "double_site", "spread7", "pair_x"]
    rows = [g(k, 500, 0.5, w) for k in ks]
    ax[2].errorbar([r["kappa_pred"] for r in rows], [r["kappa"] for r in rows], [r["kappa_sem"] for r in rows],
                   fmt=mk, label=f"omega = {w}")
    for r, k in zip(rows, ks):
        ax[2].annotate(k, (r["kappa_pred"], r["kappa"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
ax[2].plot([0, 7], [0, 7], "k:", lw=.8)
ax[2].set_xlabel("mean-field prediction")
ax[2].set_ylabel("measured kappa per unit")
ax[2].set_title("Acceleration response, m0 = 500\n(kappa > 0: pushes along the acceleration)")
ax[2].legend(fontsize=7)
fig.suptitle("EX-GC17 (EXPLORATORY): drag and acceleration response of recoil consumers — our construction", fontsize=10)
fig.tight_layout()
out = os.path.join(ROOT, "paper", "figures", "exgc17_drag.png")
fig.savefig(out, dpi=130)
print(out)
