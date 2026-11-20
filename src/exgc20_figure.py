"""EX-GC20/21 figure: kicked-particle speed for three G media and no-push control; fall in the wave medium; knot lag."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def L(d, c):
    return json.load(open(os.path.join(ROOT, "results", d, c + ".json")))["rows"]


fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for c, lab, col in (("kick_tauJ0.25", "diffusive G (D 1)", "#1f77b4"), ("kick_tauJ2", "tau_J 2 (D_eff 8)", "#ff7f0e"),
                    ("kick_tauJ20", "tau_J 20 (wave-like, D_eff 80)", "#d62728"),
                    ("kick_tauJ20_nopush", "wave-like, recoil push off", "#7f7f7f")):
    r = L("EX-GC20", c)
    ax[0].plot([x["t"] for x in r], [-x["Vcore"][2] for x in r], color=col, label=lab)
ax[0].set_xlabel("t")
ax[0].set_ylabel("speed after a kick of 0.05")
ax[0].set_title("Does the kicked particle coast?")
ax[0].legend(fontsize=7)
for c, lab, col in (("fall_tauJ20_K0.2_jinit", "K 0.2", "#d62728"), ("fall_tauJ20_K0.4_jinit", "K 0.4", "#ff7f0e")):
    r = L("EX-GC20", c)
    ax[1].plot([x["t"] for x in r], [-x["Vcore"][2] for x in r], color=col, label=lab + " speed")
    ax2 = ax[1]
ax[1].plot([0, 19], [0, 0.002 * 19], "k:", lw=.8, label="push/mass 0.002 t")
ax[1].set_xlabel("t")
ax[1].set_ylabel("downward speed")
ax[1].set_title("Fall in the wave-like medium")
ax[1].legend(fontsize=7)
for d, c, lab, col in (("EX-GC21", "knot_kfr5_lam1", "4x sweeps", "#2ca02c"), ("EX-GC21", "knot_kfr5_lam4", "4x sweeps, lambda 4", "#9467bd"),
                       ("EX-GC21", "knot_flowadv", "knot carried by the flow", "#8c564b")):
    r = L(d, c)
    ax[2].plot([x["t"] for x in r], [x["Kcen"][2] for x in r], color=col, label=lab + " (knot)")
r = L("EX-GC21", "knot_flowadv")
ax[2].plot([x["t"] for x in r], [x["Wcen"][2] for x in r], "k-", lw=2, alpha=.5, label="core (all runs)")
ax[2].set_xlabel("t")
ax[2].set_ylabel("z position")
ax[2].set_title("Knot following the falling core")
ax[2].legend(fontsize=7)
fig.suptitle("EX-GC20/21 (EXPLORATORY): wave-like G and faster knot — our construction", fontsize=10)
fig.tight_layout()
out = os.path.join(ROOT, "paper", "figures", "exgc20_wave.png")
fig.savefig(out, dpi=130)
print(out)
